"""Equilibrium dashboard for evolutive NICE runs (recorder plot file).

For the time selected in the dashboard it shows:

- Flux map: psi contours, the plasma boundary now (red) and at the start of
  the run (grey, dashed), the path of the magnetic axis so far (orange),
  X/O-points, wall, vessel and coils, with equal R and Z scales.
- Global quantities over the whole run, one quantity per plot, each axis
  fitted to its own quantity, the value at the selected time in the title and
  a cursor at that time: |Ip|, psi at the boundary, loop voltage, li_3, R/Z
  of the geometric centre and of the magnetic axis, shape, beta_pol, beta_N,
  q_95, q_axis. Steps where NICE flagged its output as invalid
  (code.output_flag = -1) are marked with red crosses.
- Coil currents (one plot per coil, plus the VS1 converter current) with
  their limits, every coil as a fraction of its limit at the selected time,
  and coil voltages.
- p', ff', q and j_phi against normalised psi, now and at the start.
- A status line with the key numbers and warnings.

Clicking on a time plot selects that time. The plots wrap to the width of
the window.

Use it for any recorder receiving NICE's equilibrium and pf_active, e.g.

    recorder_equilibrium.config: ${PDS_REPO}/visualization/nice_evo.py
    recorder_equilibrium.md: "wall=imas:hdf5?path=... pf_active=imas:hdf5?path=..."

(``md`` is optional: it only adds wall, vessel and coils to the flux map.)

What the recorder stores (one Zarr group each):

- ``equilibrium``: a few packed arrays per time slice -- the scalars, the
  boundary resampled to N_OUTLINE points, up to N_CRIT X/O-points and the
  profiles on a fixed psi_N grid.
- ``pf_active``: coil currents and voltages.
- ``flux_map``: psi on NICE's GGD mesh at every slice.
- ``flux_grid``: the mesh itself, stored only when it changes.

Every array keeps the same size from one message to the next, so the recorder
only appends. (When a size changes -- e.g. a boundary with a different number
of points, as stored by nice_inv.py -- the recorder rewrites the whole group
from scratch, which gets slower and slower on long runs.)
"""

import logging
import re

import matplotlib.tri as mtri
import numpy as np
import panel as pn
import param
import xarray as xr
from bokeh.events import Tap
from bokeh.models import (
    BoxZoomTool,
    ColorBar,
    ColumnDataSource,
    DataRange1d,
    FactorRange,
    HoverTool,
    LinearColorMapper,
    PanTool,
    Range1d,
    ResetTool,
    SaveTool,
    Span,
    WheelZoomTool,
)
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from matplotlib.figure import Figure

from imas_muscle3.visualization.base_plotter import BasePlotter
from imas_muscle3.visualization.base_state import BaseState

# Picked up by name (not called here) by the dashboard's recorder viewer.
from imas_muscle3.visualization.md_loader import init_state  # noqa: F401

logger = logging.getLogger()

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
N_OUTLINE = 201                       # boundary points stored per slice
N_CRIT = 4                            # X-points and O-points stored per slice
PSI_N = np.linspace(0.0, 1.0, 101)    # grid of the stored profiles

# Coil current limits [A] by short coil name (the name in brackets at the end
# of pf_active.coil.name, or the whole name). Same values as the controller's
# limit list; coils not listed get no limit. Check them for your machine
# description.
I_MAX = {
    "CS3U": 45e3, "CS2U": 45e3, "CS1U": 45e3, "CS1L": 45e3, "CS2L": 45e3, "CS3L": 45e3,
    "PF1": 48e3, "PF2": 55e3, "PF3": 55e3, "PF4": 55e3, "PF5": 55e3,
}
I_VS1_MAX = 22.5e3   # VS1 converter: I_PF2 + I_PF3 - I_PF4 - I_PF5 [A]

# Packed layout of the "equilibrium" group
SCALARS = (
    "ip", "li_3", "beta_pol", "beta_tor_norm", "q_95", "q_axis", "psi_axis", "psi_boundary",
    "r_axis", "z_axis", "r_geo", "z_geo", "minor_radius", "elongation",
    "triangularity_upper", "triangularity_lower", "output_flag",
)
PROFILES = ("dpressure_dpsi", "f_df_dpsi", "q", "j_phi", "pressure")
CRITICAL = ("x_r", "x_z", "o_r", "o_z")

# Meshes already stored by this recorder process, and the last time seen
# (see State._extract_flux_map)
_GRIDS_SEEN = set()
_LAST_TIME = [None]


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def _float(node):
    """A 0D IDS value as float, NaN when empty (IMAS stores -9e40)."""
    try:
        v = float(node)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) and v > -1e39 else np.nan


def _array(node):
    """An IDS array as a flat float numpy array (empty -> size 0)."""
    try:
        a = np.asarray(node, dtype=float).ravel()
    except (TypeError, ValueError):
        return np.array([])
    return np.where(a < -1e39, np.nan, a)


def _fixed(values, n):
    """values truncated or NaN-padded to length n."""
    out = np.full(n, np.nan)
    values = np.asarray(values, dtype=float).ravel()[:n]
    out[: values.size] = values
    return out


def _resample_closed(r, z, n):
    """Closed curve resampled to n points evenly spaced along its length."""
    ok = np.isfinite(r) & np.isfinite(z)
    r, z = r[ok], z[ok]
    if r.size < 3:
        return np.full(n, np.nan), np.full(n, np.nan)
    if np.hypot(r[-1] - r[0], z[-1] - z[0]) > 1e-9:
        r, z = np.append(r, r[0]), np.append(z, z[0])
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(r), np.diff(z)))])
    s_new = np.linspace(0.0, s[-1], n)
    return np.interp(s_new, s, r), np.interp(s_new, s, z)


def _short_name(name):
    """'Central Solenoid 3U (CS3U)' -> 'CS3U'."""
    m = re.search(r"\(([^()]+)\)\s*$", str(name))
    return m.group(1).strip() if m else str(name).strip()


def _finite(*arrays):
    """Flatten the arrays and keep the points where all are finite."""
    arrays = [np.ravel(np.asarray(a, dtype=float)) for a in arrays]
    ok = np.logical_and.reduce([np.isfinite(a) for a in arrays])
    return [a[ok] for a in arrays]


def _one_per_time(obj):
    """Sorted by time, one entry per time (the last one written)."""
    index = obj.indexes["time"]
    if not index.is_unique:
        obj = obj.drop_duplicates("time", keep="last")
        index = obj.indexes["time"]
    return obj if index.is_monotonic_increasing else obj.sortby("time")


def _grid_key(r, z):
    """Cheap fingerprint of a mesh."""
    half = r.size // 2
    return (r.size, float(np.sum(r)), float(np.sum(z)), float(np.sum(r[:half])), float(np.sum(z[:half])))


# ----------------------------------------------------------------------------
# What the recorder stores
# ----------------------------------------------------------------------------
class State(BaseState):
    def extract(self, message):
        name = message.metadata.name
        if name == "equilibrium":
            self._extract_equilibrium(message)
        elif name == "pf_active":
            self._extract_pf_active(message)

    def _store(self, key, new):
        current = self.data.get(key)
        if current is None:
            self.data[key] = new
        else:
            combined = xr.concat([current, new], dim="time", join="outer")
            self.data[key] = combined.drop_duplicates("time", keep="last")

    # -- pf_active ------------------------------------------------------------
    def _extract_pf_active(self, ids):
        times = _array(ids.time)
        if times.size == 0:
            return
        names = [str(c.name.value) for c in ids.coil]
        currents = np.full((times.size, len(names)), np.nan)
        voltages = np.full((times.size, len(names)), np.nan)
        for k, coil in enumerate(ids.coil):
            currents[:, k] = _fixed(_array(coil.current.data), times.size)
            voltages[:, k] = _fixed(_array(coil.voltage.data), times.size)
        self._store(
            "pf_active",
            xr.Dataset(
                {
                    "currents": (("time", "coil"), currents),
                    "voltages": (("time", "coil"), voltages),
                },
                coords={"time": times, "coil": names},
            ),
        )

    # -- equilibrium ----------------------------------------------------------
    def _extract_equilibrium(self, ids):
        # A live stream sends one slice per message; an inverse loop can send
        # a whole trace. Handle both.
        flags = _array(ids.code.output_flag)
        times = _array(ids.time)
        for i, ts in enumerate(ids.time_slice):
            t = _float(ts.time)
            if np.isnan(t) and i < times.size:
                t = times[i]
            if np.isnan(t):
                continue
            flag = flags[i] if i < flags.size else np.nan
            self._extract_slice(ts, t, flag)
            self._extract_flux_map(ts, t)

    def _extract_slice(self, ts, t, flag):
        gq = ts.global_quantities
        b = ts.boundary
        prof = ts.profiles_1d

        psi_axis = _float(gq.psi_axis)
        psi_b = _float(b.psi)
        if np.isnan(psi_b):
            psi_b = _float(gq.psi_boundary)
        scalars = {
            "ip": _float(gq.ip),
            "li_3": _float(gq.li_3),
            "beta_pol": _float(gq.beta_pol),
            "beta_tor_norm": _float(gq.beta_tor_norm),
            "q_95": _float(gq.q_95),
            "q_axis": _float(gq.q_axis),
            "psi_axis": psi_axis,
            "psi_boundary": psi_b,
            "r_axis": _float(gq.magnetic_axis.r),
            "z_axis": _float(gq.magnetic_axis.z),
            "r_geo": _float(b.geometric_axis.r),
            "z_geo": _float(b.geometric_axis.z),
            "minor_radius": _float(b.minor_radius),
            "elongation": _float(b.elongation),
            "triangularity_upper": _float(b.triangularity_upper),
            "triangularity_lower": _float(b.triangularity_lower),
            "output_flag": float(flag),
        }

        boundary = np.array(_resample_closed(_array(b.outline.r), _array(b.outline.z), N_OUTLINE))

        crit = {k: [] for k in CRITICAL}
        for node in ts.contour_tree.node:
            kind = _float(node.critical_type)
            prefix = "x" if kind == 1 else "o" if kind in (0, 2) else None
            if prefix:
                crit[f"{prefix}_r"].append(_float(node.r))
                crit[f"{prefix}_z"].append(_float(node.z))
        critical = np.array([_fixed(crit[k], N_CRIT) for k in CRITICAL])

        # Profiles on a fixed psi_N grid
        psi = _array(prof.psi)
        profiles = np.full((len(PROFILES), PSI_N.size), np.nan)
        if psi.size >= 2:
            psi_0 = psi_axis if np.isfinite(psi_axis) else psi[0]
            psi_1 = psi_b if np.isfinite(psi_b) else psi[-1]
            if psi_1 != psi_0:
                psin = (psi - psi_0) / (psi_1 - psi_0)
                for j, key in enumerate(PROFILES):
                    values = _array(getattr(prof, key))
                    if values.size != psi.size:
                        continue
                    ok = np.isfinite(psin) & np.isfinite(values)
                    if ok.sum() >= 2:
                        order = np.argsort(psin[ok])
                        profiles[j] = np.interp(PSI_N, psin[ok][order], values[ok][order],
                                                left=np.nan, right=np.nan)

        self._store(
            "equilibrium",
            xr.Dataset(
                {
                    "scalars": (("time", "quantity"), [[scalars[k] for k in SCALARS]]),
                    "boundary": (("time", "rz", "point"), [boundary]),
                    "critical": (("time", "critical_kind", "crit"), [critical]),
                    "profiles": (("time", "profile", "psi_n"), [profiles]),
                },
                coords={
                    "time": [t],
                    "quantity": list(SCALARS),
                    "critical_kind": list(CRITICAL),
                    "profile": list(PROFILES),
                    "psi_n": PSI_N,
                },
            ),
        )

    def _extract_flux_map(self, ts, t):
        """psi on NICE's GGD mesh; the mesh itself only when it changes."""
        if not (len(ts.ggd) and len(ts.ggd[0].psi) and len(ts.ggd[0].r) and len(ts.ggd[0].z)):
            return
        g = ts.ggd[0]
        r, z, psi = _array(g.r[0].values), _array(g.z[0].values), _array(g.psi[0].values)
        if not (r.size and r.size == z.size == psi.size):
            return
        self._store("flux_map", xr.Dataset({"psi": (("time", "grid_point"), [psi])}, coords={"time": [t]}))
        # Time going back = the recorder starts a new store (next outer-loop
        # iteration), which needs its own copy of the mesh.
        if _LAST_TIME[0] is not None and t < _LAST_TIME[0]:
            _GRIDS_SEEN.clear()
        _LAST_TIME[0] = t
        key = _grid_key(r, z)
        if key not in _GRIDS_SEEN:
            _GRIDS_SEEN.add(key)
            self._store(
                "flux_grid",
                xr.Dataset(
                    {"grid_r": (("time", "grid_point"), [r]), "grid_z": (("time", "grid_point"), [z])},
                    coords={"time": [t]},
                ),
            )


# ----------------------------------------------------------------------------
# What the dashboard shows
#
# Plain Bokeh figures, updated in place: new data replaces the columns of a
# few shared data sources, a new selected time moves the cursors and redraws
# the flux map, the profiles and the plot titles. Each small plot shows one
# quantity, and its axis follows that quantity (autoscaled by the browser).
# The plots sit in boxes that wrap to the width of the window, so the layout
# works at any window size or browser zoom.
# ----------------------------------------------------------------------------
TOOLS = "pan,box_zoom,wheel_zoom,reset,save"
BLUE = "#1f77b4"

# Global quantities: (column, title, unit, smallest axis span shown)
TIME_PLOTS = (
    ("ip_ma", "|Ip|", "MA", 0.02),
    ("psi_boundary", "ψ at the boundary", "Wb", 0.02),
    ("loop_voltage", "Loop voltage -dψ_b/dt", "V", 0.2),
    ("li_3", "li_3", "", 0.005),
    ("z_geo", "Z geometric centre", "m", 0.005),
    ("z_axis", "Z magnetic axis", "m", 0.005),
    ("r_geo", "R geometric centre", "m", 0.005),
    ("r_axis", "R magnetic axis", "m", 0.005),
    ("minor_radius", "Minor radius", "m", 0.005),
    ("elongation", "Elongation", "", 0.005),
    ("triangularity_upper", "Upper triangularity", "", 0.005),
    ("triangularity_lower", "Lower triangularity", "", 0.005),
    ("beta_pol", "β_pol", "", 0.002),
    ("beta_tor_norm", "β_N", "", 0.002),
    ("q_95", "q_95", "", 0.02),
    ("q_axis", "q_axis", "", 0.02),
)
# Coil tabs: (stored variable, tab title, unit, divisor, smallest axis span)
COIL_PLOTS = (
    ("currents", "Coil currents", "kA", 1e3, 0.5),
    ("voltages", "Coil voltages", "V", 1.0, 1.0),
)
# Profiles: (name, title, unit, divisor)
PROFILE_PLOTS = (
    ("dpressure_dpsi", "p'", "Pa/Wb", 1.0),
    ("f_df_dpsi", "ff'", "T² m²/Wb", 1.0),
    ("q", "Safety factor q", "", 1.0),
    ("j_phi", "Toroidal current density", "MA/m²", 1e6),
)
PLOT_HEIGHT = 180       # small plots [px]
PLOT_MIN_WIDTH = 300    # narrowest small plot before its box wraps to a new row [px]
FLUX_MAP_HEIGHT = 640   # flux-map frame [px]; the width follows the machine


def _title(title, unit):
    return title + (f" [{unit}]" if unit else "")


def _fmt(value, unit=""):
    if value is None or not np.isfinite(value):
        return "n/a"
    return f"{value:.4g}" + (f" {unit}" if unit else "")


def _pane(fig, min_width=PLOT_MIN_WIDTH):
    """A Bokeh figure as an item of a wrapping box: at least min_width wide,
    sharing the row with the others. (The Column is needed: a Bokeh pane
    alone is replaced by the figure, which ignores the flex style.)"""
    return pn.Column(pn.pane.Bokeh(fig, sizing_mode="stretch_width"), sizing_mode="stretch_width",
                     min_width=min_width, styles={"flex": f"1 1 {min_width}px"}, margin=0)


def _box(*objects):
    return pn.FlexBox(*objects, flex_wrap="wrap", sizing_mode="stretch_width")


def _plot_box(panes, min_width=PLOT_MIN_WIDTH):
    """Wrapping box of plots. A few empty items of the same size close the
    last row, so that its plots keep the width of the others."""
    fillers = [pn.Spacer(height=0, min_width=min_width, sizing_mode="stretch_width", margin=0,
                         styles={"flex": f"1 1 {min_width}px"}) for _ in range(5)]
    return _box(*panes, *fillers)


def _small_figure(x_range):
    fig = figure(height=PLOT_HEIGHT, sizing_mode="stretch_width", tools=TOOLS, x_range=x_range,
                 y_range=DataRange1d(only_visible=True, range_padding=0.1), toolbar_location="right")
    fig.toolbar.autohide = True
    fig.toolbar.logo = None
    fig.toolbar.active_scroll = None   # the mouse wheel scrolls the page
    fig.title.text_font_size = "10pt"
    fig.min_border_left = 55
    return fig


def _span_points(cols, min_spans):
    """Two invisible points per quantity that keep its axis at least
    min_span tall, centred on the data: an almost constant quantity would
    otherwise be blown up to its numerical noise."""
    t = cols["time"]
    out = {"time": np.full(2, t[-1] if len(t) else 0.0)}
    for name, span in min_spans.items():
        v = np.asarray(cols.get(name, []), dtype=float)
        v = v[np.isfinite(v)]
        if not v.size:
            out[name] = np.zeros(2)
            continue
        mid = 0.5 * (float(v.min()) + float(v.max()))
        half = 0.5 * span if np.ptp(v) < span else 0.0
        out[name] = np.array([mid - half, mid + half])
    return out


class Plotter(BasePlotter):
    levels = param.Integer(default=20, bounds=(2, 100), doc="Number of contour levels")

    # -- data access ----------------------------------------------------------
    # The viewer reads the Zarr stores while the recorder writes them, so a
    # read can fail for one refresh: the helpers return None instead of
    # raising, and the plots keep their previous content.
    def _var(self, group, name):
        """One variable over the whole run (sorted, one entry per time)."""
        ds = self.active_state.data.get(group)
        if ds is None or name not in ds:
            return None
        try:
            da = _one_per_time(ds[name])
            return da if da.sizes.get("time", 0) else None
        except Exception:
            logger.debug("%s/%s not readable", group, name, exc_info=True)
            return None

    def _load(self, da, where="all"):
        """da loaded: all times, the selected time ("now") or the first one."""
        if da is None:
            return None
        try:
            if where == "now":
                da = da.sel(time=self.time, method="nearest")
            elif where == "first":
                da = da.isel(time=0)
            return da.load()
        except Exception:
            logger.debug("not readable at t=%s", self.time, exc_info=True)
            return None

    # -- layout ---------------------------------------------------------------
    def get_dashboard(self):
        self._sc = None     # scalars (time, quantity) as of the last data update
        self._seen_data = None
        self._spans = []    # time cursors of all the time plots
        self._x_range = DataRange1d(only_visible=True, range_padding=0.02)

        globals_box = self._build_time_plots()
        coil_tabs = [self._build_coil_tab(*spec) for spec in COIL_PLOTS]
        profiles_box = self._build_profile_plots()
        flux_map = self._build_flux_map()

        levels = pn.widgets.IntSlider.from_param(self.param.levels, name="Contour levels", width=300)
        self._status_pane = pn.pane.Markdown("**Waiting for data...**", sizing_mode="stretch_width")
        tabs = pn.Tabs(
            ("Global quantities", globals_box),
            *coil_tabs,
            ("Profiles", profiles_box),
            sizing_mode="stretch_width",
            dynamic=True,
        )
        hint = pn.pane.Markdown("*Click on a time plot to select that time (tick Live View to follow "
                                "the run again). Plot tools appear on hover; the reset tool (↻) brings "
                                "back the automatic axes after a zoom.*",
                                styles={"font-size": "small"}, margin=(0, 10))
        return pn.Column(
            self._status_pane,
            _box(
                pn.Column(levels, pn.pane.Bokeh(flux_map)),
                pn.Column(hint, tabs, sizing_mode="stretch_width", min_width=2 * PLOT_MIN_WIDTH,
                          styles={"flex": f"1 1 {2 * PLOT_MIN_WIDTH}px"}),
            ),
            sizing_mode="stretch_width",
        )

    def _add_time_tools(self, fig, line, label, tooltip):
        """Time cursor, hover readout and click-to-select-time of a time plot."""
        span = Span(location=self.time, dimension="height", line_color="black", line_dash="dotted",
                    line_width=1.5)
        fig.add_layout(span)
        self._spans.append(span)
        fig.add_tools(HoverTool(renderers=[line], mode="vline",
                                tooltips=[("t", "@time{0.0000} s"), (label, tooltip)]))
        fig.on_event(Tap, self._on_tap)

    def _build_time_plots(self):
        columns = ["time"] + [c for c, *_ in TIME_PLOTS]
        self._src_scalars = ColumnDataSource({c: [] for c in columns})
        self._src_bad = ColumnDataSource({c: [] for c in columns})
        self._src_span = ColumnDataSource({c: [] for c in columns})
        self._time_figs = {}
        for column, title, unit, _ in TIME_PLOTS:
            fig = _small_figure(self._x_range)
            fig.title.text = _title(title, unit)
            line = fig.line("time", column, source=self._src_scalars, line_width=2, color=BLUE)
            fig.scatter("time", column, source=self._src_bad, marker="x", size=9, color="red", line_width=2)
            fig.scatter("time", column, source=self._src_span, alpha=0)
            self._add_time_tools(fig, line, title, f"@{{{column}}}{{0.0000}} {unit}")
            self._time_figs[column] = (fig, title, unit)
        return _plot_box([_pane(f) for f, *_ in self._time_figs.values()])

    def _build_coil_tab(self, kind, tab_title, unit, divisor, min_span):
        """One small plot per coil; the plots are made once the coil names are
        known (first data update, see _set_coil_figures)."""
        coils = {
            "unit": unit, "divisor": divisor, "min_span": min_span, "figs": {},
            "src": ColumnDataSource({"time": []}), "span": ColumnDataSource({"time": []}),
            "box": _box(pn.pane.Markdown("**Waiting for data...**")), "first": [],
        }
        self.__dict__.setdefault("_coils", {})[kind] = coils
        if kind == "currents":
            self._src_bars = ColumnDataSource({"coil": [], "fraction": [], "color": [], "label": []})
            self._bars = figure(height=120, sizing_mode="stretch_width", y_range=FactorRange(),
                                x_range=Range1d(0, 1.1), tools="hover,save", toolbar_location=None,
                                title="|I| / I_max at the selected time (limit = red line)",
                                tooltips=[("coil", "@coil"), ("|I|/I_max", "@fraction{0.000}"), ("I", "@label")])
            self._bars.hbar(y="coil", right="fraction", height=0.7, color="color", source=self._src_bars)
            self._bars.add_layout(Span(location=1.0, dimension="height", line_color="red", line_dash="dashed",
                                       line_width=2))
            self._bars.title.text_font_size = "10pt"
            # A whole row of the box (a figure above the box, instead, breaks
            # the layout of the plots in the box)
            coils["first"] = [pn.Column(pn.pane.Bokeh(self._bars, sizing_mode="stretch_width"),
                                        sizing_mode="stretch_width", styles={"flex": "1 1 100%"}, margin=0)]
        return tab_title, coils["box"]

    def _set_coil_figures(self, kind, names):
        coils = self._coils[kind]
        if list(coils["figs"]) == names:
            return
        coils["figs"] = {}
        for name in names:
            fig = _small_figure(self._x_range)
            fig.title.text = f"{name} [{coils['unit']}]"
            line = fig.line("time", name, source=coils["src"], line_width=2,
                            color="black" if name == "VS1" else BLUE)
            fig.scatter("time", name, source=coils["span"], alpha=0)
            if kind == "currents":
                limit = (I_VS1_MAX if name == "VS1" else I_MAX.get(name, np.nan)) / coils["divisor"]
                if np.isfinite(limit):
                    # Not part of the autoscale: visible when the current gets near it
                    for sign in (1, -1):
                        fig.add_layout(Span(location=sign * limit, dimension="width", line_color="red",
                                            line_dash="dashed", line_width=1.5))
            self._add_time_tools(fig, line, name, f"@{{{name}}}{{0.000}} {coils['unit']}")
            coils["figs"][name] = fig
        coils["box"].objects = coils["first"] + _plot_box([_pane(f) for f in coils["figs"].values()]).objects

    def _build_profile_plots(self):
        columns = ["psi_n"] + [f"{p}_{w}" for p, *_ in PROFILE_PLOTS for w in ("now", "first")]
        self._src_profiles = ColumnDataSource({c: [] for c in columns})
        self._profile_figs = []
        for name, title, unit, _ in PROFILE_PLOTS:
            fig = figure(height=280, sizing_mode="stretch_width", tools=TOOLS, toolbar_location="right",
                         title=_title(title, unit))
            fig.toolbar.autohide = True
            fig.toolbar.logo = None
            fig.toolbar.active_scroll = None
            fig.title.text_font_size = "10pt"
            fig.line("psi_n", f"{name}_first", source=self._src_profiles, line_width=2, color="grey",
                     line_dash="dashed", legend_label="start of the run")
            line = fig.line("psi_n", f"{name}_now", source=self._src_profiles, line_width=2, color=BLUE,
                            legend_label="selected time")
            fig.add_tools(HoverTool(renderers=[line], mode="vline",
                                    tooltips=[("ψ_N", "@psi_n{0.00}"), (title, f"@{{{name}_now}}{{0.000}}")]))
            fig.xaxis.axis_label = "ψ_N"
            legend = fig.legend[0]
            legend.orientation = "horizontal"
            legend.label_text_font_size = "8pt"
            legend.padding = 0
            legend.spacing = 12
            legend.border_line_alpha = 0
            fig.add_layout(legend, "below")
            self._profile_figs.append(fig)
        return _plot_box([_pane(f, 380) for f in self._profile_figs], 380)

    def _build_flux_map(self):
        md = getattr(self.active_state, "md", None) or {}
        limiter, vessel, quads, coil_paths = self._machine(md)

        # Window: the machine with a margin, or a default ITER-sized box; the
        # frame has the shape of the window, so R and Z have the same scale.
        r_pts = [r for r, _ in limiter + vessel + coil_paths] + [np.array([q[0], q[1]]) for q in quads]
        z_pts = [z for _, z in limiter + vessel + coil_paths] + [np.array([q[2], q[3]]) for q in quads]
        r_all, z_all = _finite(np.concatenate(r_pts), np.concatenate(z_pts)) if r_pts else ([], [])
        if len(r_all):
            x0, x1 = max(0.0, float(np.min(r_all)) - 0.5), float(np.max(r_all)) + 0.5
            y0, y1 = float(np.min(z_all)) - 0.5, float(np.max(z_all)) + 0.5
        else:
            x0, x1, y0, y1 = 1.0, 13.0, -9.0, 9.0
        frame_width = round(FLUX_MAP_HEIGHT * (x1 - x0) / (y1 - y0))

        fig = figure(frame_height=FLUX_MAP_HEIGHT, frame_width=frame_width, x_range=Range1d(x0, x1),
                     y_range=Range1d(y0, y1), title="Poloidal flux", toolbar_location="above",
                     tools=[PanTool(), WheelZoomTool(), BoxZoomTool(match_aspect=True), ResetTool(), SaveTool()])
        fig.toolbar.logo = None
        fig.xaxis.axis_label = "R [m]"
        fig.yaxis.axis_label = "Z [m]"

        self._src_contours = ColumnDataSource({"xs": [], "ys": [], "color": [], "psi": []})
        contours = fig.multi_line("xs", "ys", line_color="color", line_width=1, source=self._src_contours)
        fig.add_tools(HoverTool(renderers=[contours], tooltips=[("ψ", "@psi{0.000} Wb")]))
        self._color_mapper = LinearColorMapper(palette=Viridis256, low=0, high=1)
        fig.add_layout(ColorBar(color_mapper=self._color_mapper, title="ψ [Wb]", width=10), "right")

        if limiter:
            fig.multi_line([r for r, _ in limiter], [z for _, z in limiter], line_color="black", line_width=2)
        if vessel:
            fig.multi_line([r for r, _ in vessel], [z for _, z in vessel], line_color="grey", line_width=1)
        if quads:
            src = ColumnDataSource({k: [q[i] for q in quads] for i, k in
                                    enumerate(("left", "right", "bottom", "top", "name"))})
            coils = fig.quad(left="left", right="right", bottom="bottom", top="top", source=src,
                             fill_color="lightgrey", fill_alpha=0.5, line_color="black")
            fig.add_tools(HoverTool(renderers=[coils], tooltips=[("coil", "@name")]))
        if coil_paths:
            fig.multi_line([r for r, _ in coil_paths], [z for _, z in coil_paths], line_color="black")

        self._src_b0, self._src_axis, self._src_b, self._src_x, self._src_o = (
            ColumnDataSource({"x": [], "y": []}) for _ in range(5))
        fig.line("x", "y", source=self._src_b, color="red", line_width=3, legend_label="boundary")
        fig.line("x", "y", source=self._src_b0, color="grey", line_dash="dashed", line_width=2,
                 legend_label="boundary at start")
        fig.line("x", "y", source=self._src_axis, color="orange", line_width=2, legend_label="axis path")
        fig.scatter("x", "y", source=self._src_x, marker="x", size=12, color="black", line_width=2,
                    legend_label="X-point")
        fig.scatter("x", "y", source=self._src_o, marker="circle", size=7, color="black", legend_label="O-point")
        legend = fig.legend[0]
        legend.orientation = "horizontal"
        legend.ncols = 3
        legend.label_text_font_size = "8pt"
        legend.padding = 2
        legend.spacing = 4
        legend.border_line_alpha = 0
        fig.add_layout(legend, "below")
        return fig

    @staticmethod
    def _machine(md):
        """Limiter, vessel, coil rectangles and coil outlines from the
        machine description (empty lists without one)."""
        limiter, vessel, quads, coil_paths = [], [], [], []
        wall = md.get("wall")
        if wall is not None and len(wall.description_2d):
            d2 = wall.description_2d[0]
            for unit in d2.limiter.unit:
                limiter.append((_array(unit.outline.r), _array(unit.outline.z)))
            for unit in d2.vessel.unit:
                vessel.append((_array(unit.annular.centreline.r), _array(unit.annular.centreline.z)))
        pf = md.get("pf_active")
        if pf is not None:
            for coil in pf.coil:
                name = _short_name(coil.name.value)
                for element in coil.element:
                    geo = element.geometry
                    rect, outline, ann = geo.rectangle, geo.outline, geo.annulus
                    if _float(rect.width) > 0 and _float(rect.height) > 0:
                        r0, z0, w, h = _float(rect.r), _float(rect.z), _float(rect.width), _float(rect.height)
                        quads.append((r0 - w / 2, r0 + w / 2, z0 - h / 2, z0 + h / 2, name))
                    elif _array(outline.r).size:
                        coil_paths.append((_array(outline.r), _array(outline.z)))
                    elif _float(ann.radius_outer) > 0:
                        phi = np.linspace(0, 2 * np.pi, 17)
                        coil_paths.append((_float(ann.r) + _float(ann.radius_outer) * np.cos(phi),
                                           _float(ann.z) + _float(ann.radius_outer) * np.sin(phi)))
        return limiter, vessel, quads, coil_paths

    # -- interaction ----------------------------------------------------------
    def _on_tap(self, event):
        """Click on a time plot: select the nearest recorded time (and leave
        live view, otherwise the next data update moves the time again)."""
        times = np.asarray(self.time_slider_widget.options, dtype=float)
        if event.x is None or not times.size:
            return
        self._live_view = False
        self.time = float(times[int(np.argmin(np.abs(times - event.x)))])

    # -- new data -------------------------------------------------------------
    @param.depends("_state.data", watch=True)
    def _update_on_new_data(self):
        """New data from the recorder: refresh the time traces, then let the
        base class move the time (live view) and redraw the selected time."""
        # The recorder viewer sets the data and then triggers it again: skip
        # the second, identical, update (same datasets).
        data = self._state.data
        seen = self._seen_data
        if seen is not None and seen.keys() == data.keys() and all(seen[k] is data[k] for k in data):
            return
        self._seen_data = dict(data)   # keeps them alive, so the check above can trust "is"

        self._check_contour_cache()
        self._refresh_scalars()
        for kind, *_ in COIL_PLOTS:
            self._refresh_coils(kind)
        before = self.time
        super()._update_on_new_data()
        if self.time == before:     # otherwise _on_time already ran
            self._on_time()

    def _refresh_scalars(self):
        sc = self._load(self._var("equilibrium", "scalars"))
        if sc is None:
            return
        self._sc = sc
        t = np.asarray(sc.time.values, dtype=float)

        def q(name):
            return np.asarray(sc.sel(quantity=name).values, dtype=float)

        cols = {"time": t}
        for column, *_ in TIME_PLOTS:
            if column == "ip_ma":
                cols[column] = np.abs(q("ip")) / 1e6
            elif column == "loop_voltage":
                cols[column] = -np.gradient(q("psi_boundary"), t) if t.size > 1 else np.full(t.size, np.nan)
            else:
                cols[column] = q(column)
        self._src_scalars.data = cols
        bad = q("output_flag") == -1
        self._src_bad.data = {k: v[bad] for k, v in cols.items()}
        self._src_span.data = _span_points(cols, {c: s for c, _, _, s in TIME_PLOTS})

    def _refresh_coils(self, kind):
        da = self._load(self._var("pf_active", kind))
        if da is None:
            return
        coils = self._coils[kind]
        names = [_short_name(n) for n in da.coil.values]
        values = np.asarray(da.values, dtype=float) / coils["divisor"]
        cols = {"time": np.asarray(da.time.values, dtype=float)}
        cols.update({n: values[:, k] for k, n in enumerate(names)})
        if kind == "currents":
            vs1 = self._vs1(names, values)
            if vs1 is not None:
                cols["VS1"] = vs1
        names = [c for c in cols if c != "time"]
        coils["src"].data = cols
        coils["span"].data = _span_points(cols, dict.fromkeys(names, coils["min_span"]))
        self._set_coil_figures(kind, names)

    @staticmethod
    def _vs1(names, values):
        idx = {n: k for k, n in enumerate(names)}
        if all(c in idx for c in ("PF2", "PF3", "PF4", "PF5")):
            return values[:, idx["PF2"]] + values[:, idx["PF3"]] - values[:, idx["PF4"]] - values[:, idx["PF5"]]
        return None

    # -- selected time --------------------------------------------------------
    @param.depends("time", watch=True)
    def _on_time(self):
        for span in self._spans:
            span.location = self.time
        self._update_titles()
        for kind, *_ in COIL_PLOTS:
            self._update_coils_at_time(kind)
        self._update_flux_map()
        self._update_profiles()
        self._status_pane.object = self._status()

    @param.depends("levels", watch=True)
    def _on_levels(self):
        self._update_contours()

    def _index(self, times):
        return int(np.argmin(np.abs(np.asarray(times, dtype=float) - self.time))) if len(times) else None

    def _update_titles(self):
        data = self._src_scalars.data
        k = self._index(data["time"])
        for column, (fig, title, unit) in self._time_figs.items():
            value = data[column][k] if k is not None else None
            fig.title.text = f"{title} = {_fmt(value, unit)}"

    def _update_coils_at_time(self, kind):
        coils = self._coils[kind]
        data = coils["src"].data
        k = self._index(data["time"])
        if k is None:
            return
        bars = {"coil": [], "fraction": [], "color": [], "label": []}
        for name, fig in coils["figs"].items():
            value = float(data[name][k])
            text = f"{name} = {_fmt(value, coils['unit'])}"
            if kind == "currents":
                limit = (I_VS1_MAX if name == "VS1" else I_MAX.get(name, np.nan)) / coils["divisor"]
                if np.isfinite(limit) and np.isfinite(value):
                    frac = abs(value) / limit
                    text += f"  ({100 * frac:.0f}% of the limit)"
                    bars["coil"].append(name)
                    bars["fraction"].append(frac)
                    bars["color"].append("#d62728" if frac > 1 else "#ff7f0e" if frac > 0.9 else "#2ca02c")
                    bars["label"].append(_fmt(value, coils["unit"]))
            fig.title.text = text
        if kind == "currents":
            self._bars.y_range.factors = bars["coil"][::-1]
            self._bars.x_range.end = max([1.1] + [1.05 * f for f in bars["fraction"]])
            self._bars.height = 50 + 18 * len(bars["coil"])
            self._src_bars.data = {key: v[::-1] for key, v in bars.items()}

    # -- flux map -------------------------------------------------------------
    def _check_contour_cache(self):
        """Keep cached contours while the run only grows (live view); drop
        them when data of times already seen may have changed (another
        occurrence selected, or a store rebuilt by the recorder)."""
        ds = self._state.data.get("flux_map") if self._state is not None else None
        times = tuple(np.asarray(ds.time.values).tolist()) if ds is not None else ()
        seen = self.__dict__.get("_cache_times", ())
        if not (len(times) > len(seen) and times[: len(seen)] == seen):
            self.__dict__.setdefault("_contour_cache", {}).clear()
            self.__dict__.setdefault("_grid_cache", {}).clear()
        self.__dict__["_cache_times"] = times

    def _update_flux_map(self):
        self._update_contours()
        boundary = self._var("equilibrium", "boundary")
        first, now = self._load(boundary, "first"), self._load(boundary, "now")
        for src, b in ((self._src_b0, first), (self._src_b, now)):
            x, y = _finite(b.values[0], b.values[1]) if b is not None else ([], [])
            src.data = {"x": x, "y": y}
        if self._sc is not None:
            past = self._sc.sel(time=slice(None, self.time))
            x, y = _finite(past.sel(quantity="r_axis").values, past.sel(quantity="z_axis").values)
            self._src_axis.data = {"x": x, "y": y}
        crit = self._load(self._var("equilibrium", "critical"), "now")
        if crit is not None:
            c = {k: crit.sel(critical_kind=k).values for k in CRITICAL}
            for src, (kr, kz) in ((self._src_x, ("x_r", "x_z")), (self._src_o, ("o_r", "o_z"))):
                x, y = _finite(c[kr], c[kz])
                src.data = {"x": x, "y": y}

    def _mesh(self, t):
        """(valid-point mask, triangulation) of the mesh in use at time t."""
        grids = self._var("flux_grid", "grid_r")
        if grids is None:
            return None
        times = grids.time.values
        k = min(max(int(np.searchsorted(times, t, side="right")) - 1, 0), times.size - 1)
        cache = self.__dict__.setdefault("_grid_cache", {})
        if times[k] not in cache:
            try:
                r = np.asarray(grids.isel(time=k).load().values, dtype=float)
                z = np.asarray(self._var("flux_grid", "grid_z").isel(time=k).load().values, dtype=float)
                ok = np.isfinite(r) & np.isfinite(z)
                cache[times[k]] = (ok, mtri.Triangulation(r[ok], z[ok]))
            except Exception:
                logger.warning("Mesh at t=%s not usable", times[k], exc_info=True)
                return None
        return cache[times[k]]

    def _update_contours(self):
        psi = self._load(self._var("flux_map", "psi"), "now")
        if psi is None:
            return
        t = float(psi.time)
        cache = self.__dict__.setdefault("_contour_cache", {})
        key = (t, self.levels)
        if key not in cache:
            if len(cache) > 300:
                cache.clear()
            cache[key] = self._calc_contours(t, np.asarray(psi.values, dtype=float), self.levels)
        xs, ys, levels = cache[key]
        colors = []
        if levels:
            low, high = min(levels), max(levels)
            self._color_mapper.update(low=low, high=high if high > low else low + 1.0)
            colors = [Viridis256[round(255 * (lev - low) / ((high - low) or 1.0))] for lev in levels]
        self._src_contours.data = {"xs": xs, "ys": ys, "color": colors, "psi": levels}

    def _calc_contours(self, t, psi, levels):
        """Contour lines of psi: (list of R arrays, list of Z arrays, levels)."""
        mesh = self._mesh(t)
        if mesh is None:
            return [], [], []
        ok, tri = mesh
        if psi.size != ok.size:
            return [], [], []
        psi = psi[ok]
        if not np.isfinite(psi).all():
            # NaN breaks tricontour: drop the triangles that touch one
            bad = ~np.isfinite(psi)
            tri = mtri.Triangulation(tri.x, tri.y, tri.triangles, mask=bad[tri.triangles].any(axis=1))
            psi = np.where(bad, np.nanmean(psi), psi)
        try:
            # A private Figure, not pyplot: pyplot keeps every call in one
            # global figure that grows in a long-running dashboard.
            cs = Figure().add_subplot().tricontour(tri, psi, levels=levels)
        except (RuntimeError, ValueError):
            logger.warning("Contours skipped at t=%s", t, exc_info=True)
            return [], [], []
        xs, ys, lv = [], [], []
        for level, segs in zip(cs.levels, cs.allsegs, strict=False):
            for seg in segs:
                if len(seg) > 1:
                    xs.append(seg[:, 0])
                    ys.append(seg[:, 1])
                    lv.append(float(level))
        return xs, ys, lv

    # -- profiles -------------------------------------------------------------
    def _update_profiles(self):
        prof = self._var("equilibrium", "profiles")
        now, first = self._load(prof, "now"), self._load(prof, "first")
        if now is None:
            return
        cols = {"psi_n": PSI_N}
        for name, _, _, divisor in PROFILE_PLOTS:
            cols[f"{name}_now"] = now.sel(profile=name).values / divisor
            cols[f"{name}_first"] = (first.sel(profile=name).values / divisor if first is not None
                                     else np.full(PSI_N.size, np.nan))
        self._src_profiles.data = cols
        for fig, (_, title, unit, _) in zip(self._profile_figs, PROFILE_PLOTS, strict=True):
            fig.title.text = f"{_title(title, unit)} at t = {float(now.time):.4f} s"
            if first is not None:
                fig.legend.items[0].label = {"value": f"start ({float(first.time):.4f} s)"}

    # -- status line ------------------------------------------------------------
    def _status(self):
        sc = self._sc
        if sc is None:
            return "**Waiting for data...**"
        now = sc.sel(time=self.time, method="nearest")
        first = sc.isel(time=0)

        def v(ds, k):
            return float(ds.sel(quantity=k))

        ip_now, ip_0 = abs(v(now, "ip")), abs(v(first, "ip"))
        step = int(np.sum(sc.time.values <= float(now.time)))
        parts = [
            f"**t = {float(now.time):.4f} s** (step {step} of {sc.time.size})",
            f"|Ip| {ip_now / 1e6:.3f} MA (Δ {(ip_now - ip_0) / 1e6:+.3f})",
            f"Z geo {v(now, 'z_geo'):.3f} m (Δ {v(now, 'z_geo') - v(first, 'z_geo'):+.3f})",
            f"ψ_b {v(now, 'psi_boundary'):.2f} Wb (Δ {v(now, 'psi_boundary') - v(first, 'psi_boundary'):+.2f})",
            f"li₃ {v(now, 'li_3'):.3f}",
            f"q₉₅ {v(now, 'q_95'):.2f}",
        ]
        warnings = []
        if v(now, "output_flag") == -1:
            warnings.append("NICE output **invalid** at this step")
        n_bad = int(np.sum(sc.sel(quantity="output_flag", time=slice(None, self.time)).values == -1))
        if n_bad:
            warnings.append(f"{n_bad} invalid NICE step(s) up to this time")
        bars = self._src_bars.data
        if len(bars["fraction"]):
            k = int(np.argmax(bars["fraction"]))
            parts.append(f"max |I|/I_max {bars['fraction'][k]:.2f} ({bars['coil'][k]})")
            if bars["fraction"][k] > 1:
                warnings.append(f"**{bars['coil'][k]}** beyond its current limit")
        text = " · ".join(parts)
        if warnings:
            text += "  \n⚠ " + " · ".join(warnings)
        return text