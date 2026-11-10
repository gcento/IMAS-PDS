"""
MUSCLE3 actor between NICE and the magnetic controller. It adds two diagnostics to every
equilibrium it forwards:
- the plasma-wall gaps along the gap lines
  (time_slice[:].boundary.gap);
- the vertical velocity of the current centre, by backward difference
  (time_slice[:].global_quantities.current_centre.velocity_z).
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

def psi_on_mesh(ts) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """r, z of NICE's mesh nodes and psi on them, None if the slice has no mesh data."""
    if len(ts.ggd) == 0:
        return None
    ggd = ts.ggd[0]
    if len(ggd.r) == 0 or len(ggd.z) == 0 or len(ggd.psi) == 0:
        return None
    return np.asarray(ggd.r[0].values), np.asarray(ggd.z[0].values), np.asarray(ggd.psi[0].values)


class GapLines:

    def __init__(self) -> None:
        gaps = json.loads(GAP_FILE.read_text())["gaps"]
        self.names = [g["name"] for g in gaps]
        self.r_wall = np.array([g["r"] for g in gaps])  # start of each line, on the wall [m]
        self.z_wall = np.array([g["z"] for g in gaps])
        self.theta = np.radians([g["theta_deg"] for g in gaps])  # counter-clockwise from +R
        self.distance = np.linspace(0.0, GAP_LENGTH, N_SAMPLES)
        self.r_samples = self.r_wall[:, None] + np.cos(self.theta)[:, None] * self.distance
        self.z_samples = self.z_wall[:, None] + np.sin(self.theta)[:, None] * self.distance
        self.mesh = None
        self.triangulation = None

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
        """Gaps into ts.boundary.gap, value left empty where not found.
        The DD angle is clockwise from grad R, theta counter-clockwise."""
        ts.boundary.gap.resize(len(gaps))
        for k, gap in enumerate(ts.boundary.gap):
            gap.name = f"GAP{k + 1:02d}"
            gap.r = self.r_wall[k]
            gap.z = self.z_wall[k]
            gap.angle = -self.theta[k]
            if np.isfinite(gaps[k]):
                gap.value = gaps[k]


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


def add_diagnostics(data: bytes, gap_lines: GapLines, zdot: BackwardDifference) -> bytes:
    """Serialized equilibrium -> the same equilibrium with gaps and current-centre velocity."""
    equilibrium = IDSFactory().new("equilibrium")
    equilibrium.deserialize(data)
    output_flag = equilibrium.code.output_flag
    for i, ts in enumerate(equilibrium.time_slice):
        mesh = psi_on_mesh(ts)
        psi_b = float(ts.boundary.psi)
        if mesh is not None and abs(psi_b) < 1e30:   # NaN and empty (-9e40) fail this test
            r, z, psi = mesh
            gap_lines.write(ts, gap_lines.gaps(r, z, psi - psi_b))

        # velocity_z is written whenever the current centre is: 0 at the first step, and the
        # last one at a step NICE rejected (output_flag < 0), which repeats the previous state
        z_centre = float(ts.global_quantities.current_centre.z)
        if abs(z_centre) < 1e30:
            rejected = i < len(output_flag) and output_flag[i] < 0
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
    logger.info("diagnostics: %d gap lines from %s", len(gap_lines.r_wall), GAP_FILE.name)

    while instance.reuse_instance():
        msg = instance.receive("equilibrium_in")
        pf_active = instance.receive("pf_active_in") if pass_pf_active else None
        data = add_diagnostics(msg.data, gap_lines, zdot)
        instance.send("equilibrium_out", Message(msg.timestamp, msg.next_timestamp, data))
        if pf_active is not None:
            instance.send("pf_active_out", pf_active)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    main()
