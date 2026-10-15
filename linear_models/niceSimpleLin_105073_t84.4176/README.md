# Linear model of pulse 105073 at t = 84.417598 s (delivered under the label 105078)

Linearised plasma + circuit model delivered by the model provider (NICE linearisation),
used by this repository as an external linear model in place of the circuit model.

**Pulse correction (2026-10-14, confirmed by the user).** The provider's e-mail and the
original directory name said pulse 105078 (recorded in `provider_note` of the descriptor);
the data match pulse **105073** at t = 84.4176 s: the DINA 105073
equilibrium slice at 84.4176 s has ip = -15.0632 MA (the provider's 1.50632e7 A to the last
digit, opposite sign) and magnetic axis (6.2996, 0.6115) m (provider: (6.29704, 0.596501) m;
the 79.2176 s slice has (6.2977, 0.5958) m); the neighbouring slices have -15.0732 MA
(79.2176 s) and -15.0547 MA (89.7176 s). DINA 105078 has Ip = -7.48 MA at that time. Entry
`imas:hdf5?path=/work/projects/pds/pds-scenarios/105073/data/in`, DD 4.1.1. The descriptor's
`shot` is 105073. The user confirmed pulse 105073 on 2026-10-14; the directory was then renamed
from `linearModelSimple_105078_t84.417598` to `linearModelSimple_105073_t84.417598` (matrix
files and `linModOrder.txt` unchanged).

## Files

| File | Content |
|---|---|
| `linModSimpleA.npy`, `.txt` | A (17 x 17), continuous-time state matrix, dI/dt = A I + B U |
| `linModSimpleB.npy`, `.txt` | B (17 x 14), input matrix (14 generator voltages) |
| `linModSimpleC.npy`, `.txt` | C (9 x 17), output matrix, Y = C I |
| `linModOrder.txt` | output order: rmag, zmag, rX, zX, Ip, gaps (n_gaps = 3), Dist2ndSep |
| `linear_model.json` | descriptor read by `imas_model/configure_linear_model.m` (states, inputs, input-to-coil map, coil-to-state mapping, outputs, sign convention, operating point, notes) |

The `.txt` files hold the same matrices with 6 significant digits; the `.npy` files are
the ones used. The matrix files are untouched (provider's frame); the convention block of
the descriptor tells the wrapper how to read them.

## Sign convention (derived, 2026-10-14)

The DINA 105073 `pf_active` currents at 84.4376 s (closest sample) are CS3U -22500, CS2U
9300, CS1U = CS1L 34771, CS2L 5428, CS3L -36059, PF1 -22192, PF2 29143, PF3 20117, PF4
40774, PF5 9766, PF6 -19143, VS3U -204, VS3L 204 A; the provider's operating point (CS3U
23248.2, CS2U -9514.24, CS1 -34331.9, CS2L -6293.26, CS3L 37480.8, PF1 21873.6, PF2
-29058.6, PF3 -19884.8, PF4 -40412.3, PF5 -10061.6, PF6 19263.3, VS3 203.564, -VS3
-203.564 A) has every coil current with the opposite sign and the same magnitude within
about 3 % (CS2L 16 %). Conclusion: the matrices and the provider's operating point are
written in a frame where all currents (coils, passive, Ip) have the opposite sign to the
IMAS data (a global negation of the current sense, as between a +Ip frame and the IMAS
COCOS 17 data); geometric quantities (positions, gaps, distances) are unaffected.

Descriptor block `convention`: `state_sign -1`, `input_sign -1`, `output_sign` Ip -1 and
+1 for rmag, zmag, rX, zX, gap1..3, Dist2ndSep. Meaning: `x_file = state_sign * x`,
`u_file = input_sign * u`, `y_file(k) = output_sign(k) * y(k)` with x, u, y in the
wrapper's (IMAS) convention; `configure_linear_model` applies A unchanged,
`B_conv = (input_sign/state_sign) B` and `C_conv(k,:) = (state_sign/output_sign(k)) C(k,:)`
before the input map and the discretisation. Effect for this model: A and B unchanged, the
Ip row of C unchanged, the geometric rows of C negated. The sign of the Ip output row is
thus consistent with the IMAS data by construction. `input_sign -1` is an ASSUMPTION (the
generator voltages flip with the current sense, as power = V*I), to be confirmed; the VS3
generator-sign question below is unchanged by the frame flip (B is unchanged when both x
and u flip). The provider's second answer (2026-10-14) supports the same sense for voltages
and currents in their frame ("a positive voltage must give a positive current") and suggests
testing whether COCOS changes signs; the frame itself is still not explicitly confirmed.

## Verified content (model provider, 2026-10-14)

- States (17): the 13 coil currents `CS3U, CS2U, CS1, CS2L, CS3L, PF1 .. PF6, VS3, -VS3`
  (state 13 is the negated VS3 current, i.e. the VS3L current), then `PS_1, PS_2`
  (passive structure) and `VV_1, VV_2` (vacuum vessel).
- Inputs (14): the generator voltages `CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1, PF2, PF3,
  PF4, PF5, VS1, PF6, VS3` (VS1: star-point converter in series with PF2 .. PF5; VS3: one
  voltage for the anti-series VS3U/VS3L pair).
- I, U, Y are deviations from the equilibrium at 84.4176 s. Provider's operating point
  (provider frame): coil currents as listed above (CS1 on both CS1U and CS1L, VS3 203.564 A
  on VS3U and -203.564 A on VS3L), U = 0, rmag 6.29704 m, zmag 0.596501 m, Ip 1.50632e7 A.
  The other outputs (rX, zX, gaps, Dist2ndSep) are not used for now; their definitions
  are not needed.
- Operating point in the descriptor (wrapper/IMAS convention): the NEGATED provider coil
  currents (CS3U -23248.2, CS2U 9514.24, CS1U = CS1L 34331.9, CS2L 6293.26, CS3L -37480.8,
  PF1 -21873.6, PF2 29058.6, PF3 19884.8, PF4 40412.3, PF5 10061.6, PF6 -19263.3, VS3U
  -203.564, VS3L 203.564 A), voltages 0, rmag 6.29704 m, zmag 0.596501 m, Ip -15063200 A.
  Its origin is thereby explained: the NICE equilibrium of pulse 105073 at 84.4176 s,
  within about 3 % of the DINA 105073 coil currents (opposite frame).

## Provider's second answer and VS3 map (2026-10-14)

Model provider, 2026-10-14 (second answer; relayed by the user, translated from French):

1. On the sign frame: "the generators and the currents are linked by the circuit
   equations; normally a positive voltage must give a positive current, unless COCOS
   changes signs, which would have to be tested."
2. "VS3 and -VS3 mean VS3U and VS3L, which are currents driven by a single generator and
   placed in a circuit such that their signs are opposite."
3. On the non-mirror rows: "perhaps with the linearised model the currents end up
   decorrelated, and that is what bothers him" (the provider does not consider it normal;
   no fix provided).

Own check on the delivered B (provider frame, inputs in the confirmed generator order):
the own-generator responses B(state, its generator) are CS3U 1.522, CS2U 1.640, CS1 0.578
(per CS1U or CS1L volt), CS2L 1.631, CS3L 1.516, PF1 1.606, PF2 2.814, PF3 0.647, PF4 0.814,
PF5 0.811, PF6 0.492 A/s per V, all positive, consistent with the provider's rule. For the
VS3 pair, state `VS3` (VS3U) responds with -3192.8 A/s per V of the VS3 generator and state
`-VS3` (VS3L) with +1734.7 A/s per V. Applying the provider's rule, the positive sense of the
VS3 generator is the VS3L current direction, so the generator voltage in coil-sense
`pf_active` voltages is V(VS3L) - V(VS3U): `input_map.VS3 = {VS3U: -1, VS3L: 1}` (DERIVED,
2026-10-14; not an explicit confirmation of the loop sense). This SUPERSEDES the previous
map `{VS3U: 1, VS3L: -1}` derived from the anti-series connection and the anti-symmetric
DINA coil voltages. The two VS3 responses are also not opposite-equal (ratio -1.84), the
decorrelation the provider mentions; unresolved.

## Open questions (UNVERIFIED)

- Pulse and frame: confirm with the provider that the model is pulse 105073 at 84.4176 s
  and that its currents are in a globally negated frame relative to IMAS, including
  whether the generator voltages flip with the current sense (`input_sign -1`).
- VS3 voltage map: `input_map.VS3 = {VS3U: -1, VS3L: 1}` (generator voltage = V(VS3L) -
  V(VS3U) in coil sense) follows the provider's positive-voltage-positive-current rule
  applied to the delivered B (second answer, 2026-10-14, see above); an explicit
  confirmation of the loop sense is still pending. The previous map `{VS3U: 1, VS3L: -1}`
  (sum of the coil-sense terminal voltages around the anti-series loop, from the DINA
  anti-symmetric voltages) is superseded (the frame flip does not change B).
- VS1: no `pf_active` coil carries the VS1 generator voltage; input column 12 is not
  drivable through `pf_active` (zero row of the map; the PDS controller bridge drops the
  VS1 command anyway).
- VS3 states: A columns 12 and 13 are exact negatives, but rows 12 and 13 of A and B are
  not, so the VS3U and VS3L current deviations of a run are not mirror images. The
  provider's remark (second answer, 2026-10-14): perhaps the linearised model decorrelates
  the currents, not considered normal; no fix provided.

## Superseded history (2026-10-14)

The model was first taken as pulse 105078. The provider's Ip 1.50632e7 A was then read as
about twice the DINA 105078 value (-7.48 MA, flat between 73.6 s and 95.0 s) and treated as
a typo, and `operating_point.outputs.Ip` was set to -7478000 A by a user decision, with C
row 5 kept as delivered under the assumption of IMAS-sign Ip deviations. Both are
superseded: the pulse number was wrong, not the value (see the pulse correction above).

## How to run

See the repository `README.md`, section "External linear model (NICE linearisation)":
`run_plasmaless_imas` / `run_plasmaless_model` with `linear_model_dir` (or the environment
variable `PLASMALESS_LINEAR_MODEL_DIR`) pointing at this directory, the MUSCLE3 settings
`plasmaless.linear_model_dir` and `plasmaless.operating_point`, and the tests in
`tests/test_plasmaless_linear.m`. A PDS case matching the linearisation point must use
pulse 105073 with t0 = 84.4176 s.
