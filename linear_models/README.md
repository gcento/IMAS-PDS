# Linear plant models

Models for the `plasmaless` actor (`local_installs/plasmaless-tokamak-circuits`, branch
`IMAS_Muscle3_with_plasma`), used by `workflows/linear_plant_controller`. The model is chosen
by the setting `plant.linear_model_dir` in `cases/overrides/linear_plant_controller_<shot>.ymmsl`;
after changing it, recreate the case (`bin/pds-create-case`).

| `plant.linear_model_dir` | Model | States | Notes |
|---|---|---|---|
| `${PDS_REPO}/linear_models/plasmaless_md_dd4` | plasmaless coil+vessel circuit model, no plasma | 116 | exported by `export_plasmaless_model.m` (see below) |
| `${PDS_REPO}/linear_models/niceSimpleLin_105073_t84.4176` | NICE simplified linearisation, pulse 105073, t = 84.4176 s | 17 | copy of `linearModelSimple_105073_t84.417598` of the plasmaless repository |
| `${PDS_REPO}/linear_models/niceLin_105073_t84.4176` | NICE full linearisation, same pulse and time | 1192 | built by `make_nicelin_descriptor.py`; outputs `rX2`/`zX2` dropped; open-loop unstable |
| `""` (empty) | plasmaless circuit model built by the actor itself | 116 | the original model, from `data/md_dd4` of the plasmaless repository |

## The exported plasmaless model

`plasmaless_md_dd4` holds the same matrices as the circuit model the actor builds. Its
descriptor has `"coordinates": "absolute"`: states, inputs and outputs are absolute values,
there is no operating point, and the actor starts the model from the F_INIT coil currents at t0
(passive currents 0), exactly as with `plant.linear_model_dir: ""`. The setting
`plant.operating_point` is ignored for it.

This needs the plasmaless actor branch `IMAS_Muscle3_absolute` (fork `gcento`). The NICE models
keep the default `"coordinates": "deviation"`, so they behave as on `IMAS_Muscle3_with_plasma`.

The export script reads the folder back with the actor's own loader and checks it against the
circuit model, from the same initial coil currents: same discrete A and B, same currents over
1000 steps.

## Adding a model

A folder with:
- the matrices A, B, C (`.npy`, or whitespace-separated `.txt`): continuous time, states, inputs
  and outputs as deviations from the linearisation point;
- `linear_model.json`, in the format described in the plasmaless repository README, section
  "External linear model":
  - `states`, with `coil_states` (coil to state) and `passive_states`;
  - `inputs` (generator voltages) and `input_map` (inputs as combinations of the pf_active coil
    voltages);
  - `outputs` in {rmag, zmag, rX, zX, Ip, gap<k>, Dist2ndSep};
  - optional `operating_point` and sign `convention`.

The two NICE models use the same inputs, coil states, sign convention and operating point; a
new linearisation of the same provider can start from their descriptors.
