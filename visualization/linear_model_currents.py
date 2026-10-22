"""
Plant of linear_plant_controller (linear_model actor, muscle_linear_model_actor.m of
local_installs/plasmaless-tokamak-circuits):
- (a) CS/PF superconducting coil currents [A] (pf_active.coil(:).current)
- (b) VS resistive coil currents [kA] (coils whose name starts with "VS")
- (c) passive loop currents [kA] (pf_passive.loop(:).current)
plus the voltages carried by the same pf_active messages:
- (d) CS/PF and (e) VS (VS1, VS3) power-supply voltages [V]. The actor echoes the 14
  power-supply voltages in the coil slots (coil{k}.voltage = supply k, order SUPPLIES),
  so they are labelled with the supply names, not with the coil names of their slots
  (slots 12-14 hold VS1, PF6, VS3, the coils PF6, VS3U, VS3L). They were applied over
  the step ending at the message time: the controller's command of the previous
  exchange (the operating-point voltages, zero by default, for the first step).

Recorder ports: pf_active_in and pf_passive_in, wired as extra receivers on
plant.pf_active_o_i and plant.pf_passive_o_i (one slice per message).
"""

import logging

import holoviews as hv
import numpy as np
import panel as pn
import param
import xarray as xr
from imas_muscle3.visualization.base_plotter import BasePlotter
from imas_muscle3.visualization.base_state import BaseState

logger = logging.getLogger()

_HEIGHT = 260
_WIDTH = 900
# power supplies in the order of the voltages in pf_active (coil{k}.voltage = supply k)
SUPPLIES = ["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L", "PF1", "PF2", "PF3", "PF4",
            "PF5", "VS1", "PF6", "VS3"]


def _series(values, ntime):
    """One value per time point, NaN when the quantity is not filled."""
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size != ntime:
        return np.full(ntime, np.nan)
    return arr


def _append(current, new):
    """Concat along time; a repeated time keeps the latest values."""
    if current is None:
        return new
    combined = xr.concat([current, new], dim="time", join="outer")
    return combined.drop_duplicates("time", keep="last")


class State(BaseState):
    def extract(self, message):
        if message.metadata.name == "pf_active":
            self._extract_pf_active(message)
        elif message.metadata.name == "pf_passive":
            self._extract_pf_passive(message)

    def _extract_pf_active(self, ids):
        time = np.asarray(ids.time, dtype=float)
        if time.size == 0 or len(ids.coil) == 0:
            logger.warning("pf_active without time or coils, skipped.")
            return
        # (ncoil, ntime): one slice per message from the actor, any length works.
        current = np.array([_series(c.current.data, time.size) for c in ids.coil])
        voltage = np.array([_series(c.voltage.data, time.size) for c in ids.coil])
        names = np.array([str(c.name) for c in ids.coil])
        # voltages by power supply (by coil only if the slots are not the 14 supplies)
        supplies = SUPPLIES if len(ids.coil) == len(SUPPLIES) else names
        new = xr.Dataset(
            {
                "current": (("time", "coil"), current.T),
                "voltage": (("time", "supply"), voltage.T),
            },
            coords={"time": time, "coil": names, "supply": supplies},
        )
        self.data["pf_active"] = _append(self.data.get("pf_active"), new)

    def _extract_pf_passive(self, ids):
        time = np.asarray(ids.time, dtype=float)
        if time.size == 0 or len(ids.loop) == 0:
            logger.warning("pf_passive without time or loops, skipped.")
            return
        current = np.array([_series(lp.current, time.size) for lp in ids.loop])
        names = np.array([str(lp.name) for lp in ids.loop])
        new = xr.Dataset(
            {"current": (("time", "loop"), current.T)},
            coords={"time": time, "loop": names},
        )
        self.data["pf_passive"] = _append(self.data.get("pf_passive"), new)


class Plotter(BasePlotter):
    def get_dashboard(self):
        return pn.Column(
            hv.DynamicMap(self.plot_sc_currents),
            hv.DynamicMap(self.plot_vs_currents),
            hv.DynamicMap(self.plot_passive_currents),
            hv.DynamicMap(self.plot_sc_voltages),
            hv.DynamicMap(self.plot_vs_voltages),
        )

    def _pf_active_until_now(self, var, dim, vs):
        """(time, values, names) of var for the VS or the other entries of dim ("coil"
        or "supply"), up to the slider time; None before the first message."""
        state = self.active_state.data.get("pf_active")
        if state is None:
            return None
        names = [str(n) for n in state[dim].values]
        idx = [i for i, n in enumerate(names) if n.upper().startswith("VS") == vs]
        if not idx:
            return None
        mask = (state.time <= self.time).values
        data = state[var].isel({dim: idx, "time": mask}).transpose("time", dim)
        return state.time.values[mask], data.values, [names[i] for i in idx]

    def _overlay(self, title, var, dim, ylabel, scale, vs):
        selected = self._pf_active_until_now(var, dim, vs)
        if selected is None:
            return hv.NdOverlay(
                {"": hv.Curve(([0.0], [0.0]), "Time [s]", ylabel)}, kdims=dim
            ).opts(title="Waiting for data...", height=_HEIGHT, width=_WIDTH)
        time, values, names = selected
        curves = {
            name: hv.Curve((time, values[:, k] / scale), "Time [s]", ylabel)
            for k, name in enumerate(names)
        }
        return hv.NdOverlay(curves, kdims=dim, sort=False).opts(
            title=title,
            framewise=True,
            height=_HEIGHT,
            width=_WIDTH,
            legend_position="right",
        )

    @param.depends("time")
    def plot_sc_currents(self):
        return self._overlay(
            "(a) CS/PF superconducting coil currents",
            "current",
            "coil",
            "Current [A]",
            1.0,
            False,
        )

    @param.depends("time")
    def plot_vs_currents(self):
        return self._overlay(
            "(b) VS resistive coil currents",
            "current",
            "coil",
            "Current [kA]",
            1e3,
            True,
        )

    @param.depends("time")
    def plot_sc_voltages(self):
        return self._overlay(
            "(d) CS/PF power-supply voltages (applied over the preceding step)",
            "voltage",
            "supply",
            "Voltage [V]",
            1.0,
            False,
        )

    @param.depends("time")
    def plot_vs_voltages(self):
        return self._overlay(
            "(e) VS1 and VS3 power-supply voltages (applied over the preceding step)",
            "voltage",
            "supply",
            "Voltage [V]",
            1.0,
            True,
        )

    @param.depends("time")
    def plot_passive_currents(self):
        state = self.active_state.data.get("pf_passive")
        if state is None:
            paths, title = [], "Waiting for data..."
        else:
            mask = (state.time <= self.time).values
            time = state.time.values[mask]
            current = state.current.values[mask, :] / 1e3
            # One multi-line glyph: ~100 loops, no legend.
            paths = [
                np.column_stack([time, current[:, k]]) for k in range(current.shape[1])
            ]
            title = f"(c) Passive loop currents ({current.shape[1]} loops)"
        return hv.Path(paths, kdims=["Time [s]", "Current [kA]"]).opts(
            title=title,
            framewise=True,
            height=_HEIGHT,
            width=_WIDTH,
            color="steelblue",
            line_width=1,
        )
