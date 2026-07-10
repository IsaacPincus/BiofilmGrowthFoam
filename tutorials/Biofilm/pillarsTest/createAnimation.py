"""
Filter OpenFOAM time directories down to the ones that line up with the
overall step (e.g. 1800, 3600, ...), then compile frames saved in
./Images/ into an mp4 at a given frame rate.

Two modes, controlled by FILTER_BY_INDEX below:
  - False: assumes ./Images/ already only contains the whole-step frames
           (e.g. because you used ExtractTimeSteps in ParaView). Just
           takes every image in order.
  - True:  assumes ./Images/ contains frames for ALL timesteps, and only
           picks out the ones whose position matches the whole-step
           indices computed from the case directory.
"""

import os
import re
import glob
import imageio.v2 as imageio

STEP = 1800        # overall step size
IMAGE_DIR = "./Images/"
OUTPUT_FILE = "./B_output_filter.mp4"
FPS = 4

FILTER_BY_INDEX = True  # toggle: True = pick out only whole-step frames from a full set


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
        if abs(round(t / step) * step - t) < tol:
            entries.append((t, name))
    entries.sort(key=lambda x: x[0])
    return entries


def get_whole_step_indices(all_times, step=STEP, tol=1e-6):
    """Given a full list of timestep values, return the indices of the
    ones that are whole multiples of `step`."""
    indices = []
    for i, t in enumerate(all_times):
        if abs(round(t / step) * step - t) < tol:
            indices.append(i)
    return indices


def natural_sort_key(filename):
    """Sort filenames containing numbers in numeric order rather than
    lexicographic order (so frame_2.png comes before frame_10.png)."""
    parts = re.split(r"(\d+)", filename)
    return [int(p) if p.isdigit() else p for p in parts]


def make_video(image_dir, output_file, fps, indices=None):
    files = sorted(glob.glob(os.path.join(image_dir, "*.png")), key=natural_sort_key)
    if not files:
        raise FileNotFoundError(f"No .png files found in {image_dir}")

    if indices is not None:
        before = len(files)
        try:
            files = [files[i] for i in indices]
        except IndexError:
            raise IndexError(
                f"An index in 'indices' (max {max(indices)}) is out of range "
                f"for the {before} files found in {image_dir}. "
                "Check that the image export matches the full timestep list, "
                "not just the filtered one."
            )
        print(f"Filtered {before} frames down to {len(files)} using indices")
    else:
        print(f"Found {len(files)} frames in {image_dir}")

    with imageio.get_writer(output_file, fps=fps) as writer:
        for f in files:
            writer.append_data(imageio.imread(f))

    print(f"Wrote {len(files)} frames to {output_file} at {fps} fps")


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

    make_video(IMAGE_DIR, OUTPUT_FILE, FPS, indices=indices if FILTER_BY_INDEX else None)