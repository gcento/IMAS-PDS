# linear_plant_controller

Runs the magnetic controller `my_contr` (`controllers/my_contr`) in closed loop with a
linear plant model.

```
source ──► plant (linear_model actor) ──► magnetic_controller (my_contr)
              ▲                                                │
              └──────────── power-supply voltages ─────────────┘
```

- `source` sends the initial equilibrium and coil currents (built by `preprocess.sh`).
- `plant` integrates the chosen linear model and returns equilibrium (Ip, magnetic axis),
  coil and passive currents.
- `magnetic_controller` returns the 14 power-supply voltages.

The available models and their limits are listed in
`local_installs/plasmaless-tokamak-circuits/linear_models/README.md`.

To try a model alone before a closed-loop run, use
`local_installs/plasmaless-tokamak-circuits/utils/run_supply_sequence.m`. It gives the
open-loop response of the 12 circuit currents to a voltage on each power supply in turn,
for one or more models. To build a model for another shot or time, see
`generators/README.md` of the same repository (launchers in its `utils/`).

## Before the first run

- Run `inverse_convergence` for the same shot. Its `out_nice` gives the initial state and
  the controller references (`preprocess.source_uri`).
- The `linear_model` actor comes from `local_installs/plasmaless-tokamak-circuits`,
  branch `IMAS_Muscle3` (`setup_files/setup_plasmaless.sh`).

## Running

```bash
# create the case (runs the preprocess; MATLAB, so on a compute node)
bin/pds-create-case linear_plant_controller 105073
# run it, optionally with overlays for one-off changes
sbatch bin/pds-run-case.sbatch cases/linear_plant_controller_105073
```

## Main settings (per shot, in the case override)

| Setting | Meaning |
|---|---|
| `plant.linear_model_dir` | model folder (see the models README) |
| `plant.operating_point` | `descriptor` (model's linearisation point) or `f_init` (run's initial state) |
| `source.t_min`, `source.t_max` | time window of the initial state; the run starts at the first slice in it |
| `plant.t_end` | end of the run |
| `plant.dt`, `plant.t_interval`, `magnetic_controller.in_out_msg_dt` | exchange period, all three equal; use 0.001 with VS on |
| `magnetic_controller.KCURR` | CS/PF current controller (1/0) |
| `magnetic_controller.control_Ip` | Ip control inside KCURR; also selects KCURR's gains: 0 = plasmaless set, 1 = plasma set |
| `magnetic_controller.VS` | vertical stabilisation (needed with the unstable models) |
| `magnetic_controller.VS_control_mode` | 1 = VS1, 2 = VS3, 3 = VS1+VS3 |
| `magnetic_controller.RZIp`, `magnetic_controller.Shape_control` | RZIp and shape controllers (off by default) |
| `magnetic_controller.plasma_position` | `magnetic_axis`  or `geometric_axis` |

## Outputs (in `cases/runs/<run>/`)

- `out_plant/`: equilibrium, pf_active and pf_passive of the plant.
- `out_controller/`: pf_active with the controller's voltages.
- `simulink_log.csv`: controller internals (errors, references, VS commands, Zdot).
- Live currents for m3dash: `instances/recorder_plant/workdir/`.

## Good to know

- `my_contr` sends the 14 power-supply voltages (`coil{k}.voltage` = supply k, order
  `CS3U CS2U CS1U CS1L CS2L CS3L PF1..PF5 VS1 PF6 VS3`), the inputs of every linear
  model. Every model has the ITER circuits inside (CS1U/CS1L in series, VS3U/VS3L in
  anti-series) and returns the 14 coil currents, as NICE; `plasmaless_model` is the
  plasma-free one.
