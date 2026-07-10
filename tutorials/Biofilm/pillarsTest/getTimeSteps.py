"""
Filter OpenFOAM time directories down to the ones that line up with the
overall step (e.g. 1800, 3600, ...), dropping intermediate sub-iteration
dirs like 10800.1, 10800.2, etc.

Run standalone (lists indices from the case dir), or import
`get_whole_step_indices` inside pvpython to feed directly into
ExtractTimeSteps / reader.UpdatePipeline loops.
"""

import os

STEP = 1800  # overall step size


def get_whole_step_times(case_dir, step=STEP, tol=1e-6):
    """Return sorted list of (float time, dirname) for dirs that are
    whole multiples of `step`."""
    entries = []
    for name in os.listdir(case_dir):
        path = os.path.join(case_dir, name)
        if not os.path.isdir(path):
            continue
        try:
            t = float(name)
        except ValueError:
            continue
        # keep only multiples of step (within floating point tolerance)
        if abs(round(t / step) * step - t) < tol:
            entries.append((t, name))
    entries.sort(key=lambda x: x[0])
    return entries


def get_whole_step_indices(all_times, step=STEP, tol=1e-6):
    """Given a full list of timestep values (e.g. reader.TimestepValues
    in ParaView), return the indices of the ones that are whole
    multiples of `step`."""
    indices = []
    for i, t in enumerate(all_times):
        if abs(round(t / step) * step - t) < tol:
            indices.append(i)
    return indices


if __name__ == "__main__":
    case_dir = "./yieldStress/"  # change to your case directory if running standalone
    whole_steps = get_whole_step_times(case_dir)

    print("Whole-step time directories found:")
    for t, name in whole_steps:
        print(f"  {name}")

    all_dirs = sorted(
        (d for d in os.listdir(case_dir) if os.path.isdir(os.path.join(case_dir, d))),
        key=lambda n: (float(n) if n.replace('.', '', 1).isdigit() else -1)
    )
    all_times = []
    for d in all_dirs:
        try:
            all_times.append(float(d))
        except ValueError:
            pass
    all_times.sort()

    indices = get_whole_step_indices(all_times)
    print("\nIndices into the full sorted timestep list:")
    print(indices)