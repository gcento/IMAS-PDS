"""Build the niceLin model folder for the plasmaless actor from the NICE full linear model.

Input: the full linear model delivered for pulse 105073 at t = 84.417598 s
(linModA/B/C.npy + linModOrder.txt), and the descriptor of the simplified model of the
same provider, pulse and time (niceSimpleLin_105073_t84.4176/linear_model.json).

Output folder (read through the actor setting linear_model_dir):
- linModA.npy, linModB.npy, linModOrder.txt: copied unchanged;
- linModC_noX2.npy: C without the rows rX2, zX2. The actor (configure_linear_model.m)
  accepts only the outputs rmag, zmag, rX, zX, Ip, gap<k>, Dist2ndSep; the remaining
  outputs are those of the simplified model;
- linear_model.json: inputs, input_map, coil_states, convention, operating_point and
  units copied from the simplified model; 13 coil states (same names as the simplified
  model) followed by 1179 passive states P_1 .. P_1179.

Checks before writing: matrix sizes, and the first 13 columns of C equal to the ones of
the simplified model for every common output (same coil states, same order and sign).

Usage, once, on a compute node:
    python make_nicelin_descriptor.py SOURCE_DIR SIMPLE_DIR OUT_DIR
"""

import argparse
import datetime
import json
import shutil
from pathlib import Path

import numpy as np

N_COIL_STATES = 13
FULL_OUTPUTS = ["rmag", "zmag", "rX", "zX", "rX2", "zX2", "Ip",
                "gap1", "gap2", "gap3", "Dist2ndSep"]
DROPPED_OUTPUTS = ["rX2", "zX2"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source_dir", type=Path,
                        help="folder with linModA/B/C.npy and linModOrder.txt")
    parser.add_argument("simple_dir", type=Path,
                        help="niceSimpleLin folder (linear_model.json + matrices)")
    parser.add_argument("out_dir", type=Path, help="niceLin folder to write")
    args = parser.parse_args()

    a = np.load(args.source_dir / "linModA.npy")
    b = np.load(args.source_dir / "linModB.npy")
    c = np.load(args.source_dir / "linModC.npy")
    simple = json.loads((args.simple_dir / "linear_model.json").read_text())
    n = a.shape[0]
    m = len(simple["inputs"])
    if a.shape != (n, n) or b.shape != (n, m) or c.shape != (len(FULL_OUTPUTS), n):
        raise ValueError(f"A {a.shape}, B {b.shape}, C {c.shape} do not match {n} states, "
                         f"{m} inputs, {len(FULL_OUTPUTS)} outputs {FULL_OUTPUTS}")
    simple_states = simple["states"][:N_COIL_STATES]
    if set(simple_states) & set(simple["passive_states"]):
        raise ValueError(f"the first {N_COIL_STATES} states of the simplified model are not "
                         f"all coil states: {simple_states}")

    # same coil states: equal C coefficients on the first 13 states, output by output
    c_simple = np.load(args.simple_dir / simple["files"]["C"])
    for k, name in enumerate(simple["outputs"]):
        row_full = c[FULL_OUTPUTS.index(name), :N_COIL_STATES]
        row_simple = c_simple[k, :N_COIL_STATES]
        if not np.allclose(row_full, row_simple, rtol=1e-6, atol=0.0):
            raise ValueError(f"output {name}: the first {N_COIL_STATES} columns of C differ "
                             f"from the simplified model (max |diff| "
                             f"{np.max(np.abs(row_full - row_simple)):.3g})")
    print(f"C coil columns equal to the simplified model for {simple['outputs']}")

    outputs = [name for name in FULL_OUTPUTS if name not in DROPPED_OUTPUTS]
    if outputs != simple["outputs"]:
        raise ValueError(f"kept outputs {outputs} differ from the simplified model "
                         f"{simple['outputs']}")
    c_kept = c[[FULL_OUTPUTS.index(name) for name in outputs], :]

    passive = [f"P_{k}" for k in range(1, n - N_COIL_STATES + 1)]
    descriptor = {
        "description": ("NICE full linear model of ITER pulse 105073 at t = 84.417598 s: "
                        "dI/dt = A I + B U, Y = C I (deviations from the linearisation "
                        "point); C rows rX2, zX2 dropped"),
        "shot": simple["shot"],
        "time": simple["time"],
        "files": {"A": "linModA.npy", "B": "linModB.npy", "C": "linModC_noX2.npy"},
        "states": simple_states + passive,
        "inputs": simple["inputs"],
        "input_map": simple["input_map"],
        "coil_states": simple["coil_states"],
        "passive_states": passive,
        "outputs": outputs,
        "convention": simple["convention"],
        "operating_point": simple["operating_point"],
        "units": simple["units"],
        "provenance": (
            "matrices: NICE full linear model (directory linearModel_105073_t84.417598, "
            "linModA/B/C.npy, C rows rX2 and zX2 dropped); states 1-13, inputs, input_map, "
            "coil_states, convention, operating_point and units copied from "
            "niceSimpleLin_105073_t84.4176/linear_model.json (same provider, pulse and time; "
            "the first 13 C columns of the two models are equal): ASSUMPTION to be confirmed "
            "by the model provider; the passive states P_1 .. P_1179 are unnamed in the "
            "delivery. Generated by linear_models/make_nicelin_descriptor.py on "
            f"{datetime.date.today().isoformat()}."),
        "notes": [
            "See niceSimpleLin_105073_t84.4176/linear_model.json for the derivation of the "
            "sign frame, the VS3 input map and the operating point.",
            "Open-loop unstable (vertical mode): without vertical control the deviations "
            "diverge within a few tenths of a second.",
        ],
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("linModA.npy", "linModB.npy", "linModOrder.txt"):
        shutil.copyfile(args.source_dir / name, args.out_dir / name)
    np.save(args.out_dir / "linModC_noX2.npy", c_kept)
    (args.out_dir / "linear_model.json").write_text(json.dumps(descriptor, indent=2) + "\n")
    print(f"wrote {args.out_dir}: {n} states ({N_COIL_STATES} coil, {len(passive)} passive), "
          f"{m} inputs, {len(outputs)} outputs")


if __name__ == "__main__":
    main()
