# linear_plant_controller

`my_contr` (`controllers/my_contr`) in closed loop with a linear plant model, chosen per case:

| `plant.linear_model_dir` | Plant |
|---|---|
| `${PDS_REPO}/linear_models/plasmaless_md_dd4` | plasmaless coil+vessel circuit model, no plasma (116 states) |
| `${PDS_REPO}/linear_models/niceSimpleLin_105073_t84.4176` | NICE simplified linearisation of pulse 105073 at 84.4176 s (17 states) |
| `${PDS_REPO}/linear_models/niceLin_105073_t84.4176` | NICE full linearisation, same pulse and time (1192 states, open-loop unstable) |
| `""` | the plasmaless circuit model built by the actor itself (the original; see `linear_models/README.md`) |

The plant is always the plasmaless actor (`local_installs/plasmaless-tokamak-circuits`, branch
`IMAS_Muscle3_with_plasma`; program `plasmaless` in `workflows/lib/local_programs.ymmsl`). It
reads the model from `plant.linear_model_dir`; see `linear_models/README.md` for the models
and for adding one. The workflow is `plasmaless_controller`'s, with:
- the components renamed: `plant`, `recorder_plant`;
- the local `magnetic_controller` (`my_contr`) instead of the EasyBuild one (`KCURR_RZIp`).

## Running it

The model is set in `cases/overrides/linear_plant_controller_<shot>.ymmsl`. Keep exactly one
`plant.linear_model_dir` line active there; after changing it, recreate the case. Both commands
go through Slurm, as the preprocess and the run need compute resources:

```bash
srun -A iter --partition=all_debug --time=00:30:00 bin/pds-create-case linear_plant_controller 105073
sbatch bin/pds-run-case.sbatch cases/linear_plant_controller_105073
```

**Prerequisite.** `preprocess.sh` runs `plasmaless_controller`'s preprocess, which reads the NICE
inverse output `cases/runs/metis_nice_inverse_from_dina_<shot>/nice_out` of the same checkout.
Run that workflow for the shot first.

## Conventions

- **Exchange period.** `plant.dt = plant.t_interval = 0.002`, the period hard-coded in
  `my_contr` (`in_out_dt`).
- **Voltage order.** `env.sh` sets `MY_CONTR_VOLTAGE_ORDER=coil`, so `my_contr` sends its
  voltages in coil order, which is what the plasmaless actor reads:
  - PF6 goes in slot 12;
  - the VS3 generator voltage is split on the anti-series pair, VS3U = -VS3/2 and VS3L = +VS3/2;
  - VS1 is dropped: no coil carries it, and the `input_map` of the NICE models leaves it
    undriven.
  With NICE evolutive, `my_contr` keeps its default (power-supply order).
- **Outputs.** `out_plant` holds the plant's equilibrium, pf_active and pf_passive; for the
  NICE models the equilibrium carries the model's Ip and magnetic axis. `out_controller` holds
  the controller's voltages, and `simulink_log.csv` is written in the run directory.
- **Magnetic vs geometric axis.** The NICE models output the magnetic axis. The boundary
  geometric axis read by `my_contr` stays the reference pass-through.
