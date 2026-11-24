"""
MUSCLE3 actor between NICE and the magnetic controller. It adds three diagnostics to every
equilibrium it forwards:
- the plasma-wall gaps along the gap lines, every gap in every slice, NaN where not found
  (time_slice[:].boundary.gap);
- the vertical velocity of the current centre, by backward difference
  (time_slice[:].global_quantities.current_centre.velocity_z);
- the X-point moved from NICE's mesh node to the saddle of a local quadratic fit of psi
  (time_slice[:].contour_tree.node with critical_type 1);
- with the setting current_centre_area_weighted, the current centre recomputed area-weighted
  on the mesh (time_slice[:].global_quantities.current_centre.r/z).
"""
import logging
from pathlib import Path

import numpy as np
from imas import IDSFactory
from libmuscle import Instance, InstanceFlags, Message
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
from ymmsl.v0_2 import Operator
import json

logger = logging.getLogger()


GAP_FILE = Path(__file__).with_name("gaps.json")
GAP_LENGTH = 2.0     # [m] searched along each line, longer than any gap
N_SAMPLES = 1001     # points per segment to search on to find gap value
OUTLINE_WINDOW = (-0.1, 3.0)   # [m] along each line, where the outline crossing is searched from behind the wall point, so a gap can be negative
Z_DIVERTOR = -3.2    # [m] gaps whose wall point is lower are divertor gaps (GAP21-30, GAP59-62)

def on_mesh(ts, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """r, z of NICE's mesh nodes and `field` (psi, j_phi) on them, None if the slice has no mesh data."""
    if len(ts.ggd) == 0:
        return None
    ggd = ts.ggd[0]
    values = getattr(ggd, field)
    if len(ggd.r) == 0 or len(ggd.z) == 0 or len(values) == 0:
        return None
    return np.asarray(ggd.r[0].values), np.asarray(ggd.z[0].values), np.asarray(values[0].values)


def fit_quadratic_xpoint(ts, r: np.ndarray, z: np.ndarray, psi: np.ndarray) -> None:
    """Moves the X-point (contour_tree node with critical_type 1, a mesh node in NICE) to the
    saddle of psi ~ c0 + c1 x + c2 y + c3 x^2 + c4 x y + c5 y^2, fitted by least squares on
    the nodes within XPOINT_FIT_RADIUS. The node is left as it is when no saddle is found near it."""
    for node in ts.contour_tree.node:
        if node.critical_type == 1:
            r0, z0 = float(node.r), float(node.z)
            near = np.hypot(r - r0, z - z0) < 0.10
            if np.count_nonzero(near) < 6:
                return
            x, y = r[near] - r0, z[near] - z0
            A = np.column_stack([np.ones_like(x), x, y, x**2, x * y, y**2])
            c = np.linalg.lstsq(A, psi[near], rcond=None)[0]
            H = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])   # Hessian of the fit
            if np.linalg.det(H) >= 0:                             # not a saddle
                return
            s = -np.linalg.solve(H, c[1:3])                       # where grad psi = 0
            if np.hypot(s[0], s[1]) < 0.05:
                node.r = r0 + s[0]
                node.z = z0 + s[1]
            return


class GapLines:

    def __init__(self) -> None:
        gaps = json.loads(GAP_FILE.read_text())["gaps"]
        self.names = [g["name"] for g in gaps]
        self.r_wall = np.array([g["r"] for g in gaps])  # start of each line, on the wall [m]
        self.z_wall = np.array([g["z"] for g in gaps])
        self.theta = np.radians([g["theta_deg"] for g in gaps])  # counter-clockwise from +R
        self.divertor = self.z_wall < Z_DIVERTOR
        self.distance = np.linspace(0.0, GAP_LENGTH, N_SAMPLES)
        self.r_samples = self.r_wall[:, None] + np.cos(self.theta)[:, None] * self.distance
        self.z_samples = self.z_wall[:, None] + np.sin(self.theta)[:, None] * self.distance
        self.mesh = None
        self.triangulation = None

    def find(self, ts, mesh, psi_b: float) -> np.ndarray:
        """Gaps [m] of one slice, NaN where not found. """
        found = np.full(len(self.names), np.nan)
        r_out = np.asarray(ts.boundary.outline.r, dtype=float)
        z_out = np.asarray(ts.boundary.outline.z, dtype=float)
        if len(r_out) >= 3 and len(r_out) == len(z_out) and np.all(np.abs(np.r_[r_out, z_out]) < 1e30):
            main = ~self.divertor
            found[main] = self.outline_crossings(r_out, z_out)[main]
        if mesh is not None and abs(psi_b) < 1e30:   # NaN and empty (-9e40) fail this test
            r, z, psi = mesh
            search = self.divertor | np.isnan(found)
            found[search] = self.gaps(r, z, psi - psi_b)[search]
        return found

    def outline_crossings(self, r_out: np.ndarray, z_out: np.ndarray) -> np.ndarray:
        """Signed distance [m] from the wall point to the first crossing of the closed outline
        along each line, within OUTLINE_WINDOW, NaN where it is not crossed."""
        d_r = np.roll(r_out, -1) - r_out            
        d_z = np.roll(z_out, -1) - z_out
        u_r = np.cos(self.theta)[:, None]            # line directions
        u_z = np.sin(self.theta)[:, None]
        w_r = r_out - self.r_wall[:, None]            # wall point to segment start, lines x segments
        w_z = z_out - self.z_wall[:, None]
        den = u_r * d_z - u_z * d_r                   # 0 for a segment parallel to the line or empty
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w_r * d_z - w_z * d_r) / den         # along the line
            s = (w_r * u_z - w_z * u_r) / den         # along the segment, 0 to 1 on it
        hit = (den != 0) & (s >= 0) & (s <= 1) & (t >= OUTLINE_WINDOW[0]) & (t <= OUTLINE_WINDOW[1])
        first = np.where(hit, t, np.inf).min(axis=1)
        return np.where(np.isfinite(first), first, np.nan)

    def gaps(self, r: np.ndarray, z: np.ndarray, flux: np.ndarray) -> np.ndarray:
        """Distance [m] from the wall to the separatrix along each line, NaN where it is not crossed.
        flux is psi - psi_boundary on the mesh nodes (r, z): its sign changes on the separatrix."""
        mesh = np.column_stack([r, z])
        if self.mesh is None or not np.array_equal(mesh, self.mesh):
            self.mesh = mesh
            self.triangulation = Delaunay(mesh)     # NICE's mesh is the same at every step
        f = LinearNDInterpolator(self.triangulation, flux)(self.r_samples, self.z_samples)

        before, after = f[:, :-1], f[:, 1:]
        crosses = np.isfinite(before) & np.isfinite(after) & (np.sign(before) != np.sign(after))
        gaps = np.full(len(f), np.nan)
        for line in np.flatnonzero(crosses.any(axis=1)):
            k = np.argmax(crosses[line])             # first crossing from the wall
            d0, d1 = self.distance[k], self.distance[k + 1]
            f0, f1 = before[line, k], after[line, k]
            gaps[line] = d0 + (d1 - d0) * f0 / (f0 - f1)
        return gaps

    def write(self, ts, gaps: np.ndarray) -> None:
        """Gaps into ts.boundary.gap"""
        ts.boundary.gap.resize(len(gaps))
        for k, gap in enumerate(ts.boundary.gap):
            gap.name = f"GAP{k + 1:02d}"
            gap.r = self.r_wall[k]
            gap.z = self.z_wall[k]
            gap.angle = -self.theta[k]
            gap.value = gaps[k]


class CurrentCentre:
    """Area-weighted current centre: sum over the Delaunay triangles of the mesh of
    area * mean nodal j_phi * centroid. NICE's global_quantities.current_centre is a plain
    j_phi-weighted mean of the nodes, which depends on the mesh (16 cm between the 7340-node
    inverse and the 10776-node evolutive meshes in the limiter phase of 105073)."""

    def __init__(self) -> None:
        self.key = None

    def __call__(
        self, r: np.ndarray, z: np.ndarray, j: np.ndarray
    ) -> tuple[float, float]:
        key = (len(r), r[0], z[0], r[-1], z[-1])
        if key != self.key:  # the mesh is the same at every step
            self.key = key
            tri = self.tri = Delaunay(np.column_stack([r, z])).simplices
            self.area = 0.5 * np.abs(
                (r[tri[:, 1]] - r[tri[:, 0]]) * (z[tri[:, 2]] - z[tri[:, 0]])
                - (r[tri[:, 2]] - r[tri[:, 0]]) * (z[tri[:, 1]] - z[tri[:, 0]])
            )
            self.r_c, self.z_c = r[tri].mean(axis=1), z[tri].mean(axis=1)
        w = self.area * j[self.tri].mean(axis=1)
        return float(w @ self.r_c / w.sum()), float(w @ self.z_c / w.sum())


class BackwardDifference:
    """Time derivative of a sampled signal, (x_k - x_k-1) / (t_k - t_k-1); 0 until the
    second sample."""

    def __init__(self) -> None:
        self.t = None        # last sample
        self.x = None
        self.derivative = 0.0

    def update(self, t: float, x: float) -> float:
        """Adds the sample (t, x) and returns the derivative at t."""
        if self.t is not None and t > self.t:
            self.derivative = (x - self.x) / (t - self.t)
        self.t, self.x = t, x
        return self.derivative


def add_diagnostics(
    data: bytes,
    gap_lines: GapLines,
    zdot: BackwardDifference,
    centre: CurrentCentre | None = None,
) -> bytes:
    """Serialized equilibrium -> the same equilibrium with gaps, current-centre velocity,
    the refined X-point and, with centre, the area-weighted current centre."""
    equilibrium = IDSFactory().new("equilibrium")
    equilibrium.deserialize(data)
    output_flag = equilibrium.code.output_flag
    for i, ts in enumerate(equilibrium.time_slice):
        mesh = on_mesh(ts, "psi")
        psi_b = float(ts.boundary.psi)
        gap_lines.write(ts, gap_lines.find(ts, mesh, psi_b))
        if mesh is not None and abs(psi_b) < 1e30:   # NaN and empty (-9e40) fail this test
            r, z, psi = mesh
            fit_quadratic_xpoint(ts, r, z, psi)

        rejected = i < len(output_flag) and output_flag[i] < 0
        if centre is not None and not rejected:
            mesh_j = on_mesh(ts, "j_phi")
            if mesh_j is not None and np.isfinite(mesh_j[2]).all() and mesh_j[2].any():
                cc = ts.global_quantities.current_centre
                cc.r, cc.z = centre(*mesh_j)

        # velocity_z is written whenever the current centre is: 0 at the first step, and the
        # last one at a step NICE rejected (output_flag < 0), which repeats the previous state
        z_centre = float(ts.global_quantities.current_centre.z)
        if abs(z_centre) < 1e30:
            if not rejected:
                zdot.update(float(ts.time), z_centre)
            ts.global_quantities.current_centre.velocity_z = zdot.derivative
    return equilibrium.serialize()


def main() -> None:
    # F_INIT -> O_F: one iteration per NICE step; reuse_instance() is False once NICE closes its port
    instance = Instance(
        {
            Operator.F_INIT: ["equilibrium_in", "pf_active_in"],
            Operator.O_F: ["equilibrium_out", "pf_active_out"],
        },
        InstanceFlags.STATE_NOT_REQUIRED_FOR_NEXT_USE,   # it keeps only a cache and the last sample
    )
    pass_pf_active = instance.is_connected("pf_active_in")   # forwarded unchanged
    gap_lines = GapLines()
    zdot = BackwardDifference()
    centre = (
        CurrentCentre()
        if instance.get_setting("current_centre_area_weighted", "bool", default=False)
        else None
    )
    logger.info("diagnostics: %d gap lines from %s", len(gap_lines.r_wall), GAP_FILE.name)

    while instance.reuse_instance():
        msg = instance.receive("equilibrium_in")
        pf_active = instance.receive("pf_active_in") if pass_pf_active else None
        data = add_diagnostics(msg.data, gap_lines, zdot, centre)
        instance.send("equilibrium_out", Message(msg.timestamp, msg.next_timestamp, data))
        if pf_active is not None:
            instance.send("pf_active_out", pf_active)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    main()
