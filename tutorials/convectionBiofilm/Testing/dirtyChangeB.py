import re
import numpy as np

# ---------------- USER SETTINGS ----------------
input_file = "0/B"
output_file = "0/B"
c = 8000.0
clip_min = 0.0
clip_max = 1.0
# ----------------------------------------------

with open(input_file, "r") as f:
    text = f.read()

# Match nonuniform internalField
pattern = r"(internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\()(.*?)(\)\s*;)"
match = re.search(pattern, text, re.S)

if not match:
    raise RuntimeError("Could not find nonuniform internalField in 0/B")

header, values_block, footer = match.groups()

# Parse numbers (robust against newlines/spaces)
values = np.fromstring(values_block, sep=' ')

if values.size == 0:
    raise RuntimeError("Failed to parse field values")

# ---- STEP 1: clip original field to [0,1]
values_clipped = np.clip(values, clip_min, clip_max)

# ---- STEP 2: transform
new_values = c * (1.0 - values_clipped)

# ---- rebuild OpenFOAM list format
new_block = "\n" + "\n".join(f"{v:.10g}" for v in new_values) + "\n"

# ---- replace in file
new_text = re.sub(
    pattern,
    lambda m: m.group(1) + new_block + m.group(3),
    text,
    flags=re.S
)

with open(output_file, "w") as f:
    f.write(new_text)

print("✔ B field clipped to [0,1], transformed, and written back.")