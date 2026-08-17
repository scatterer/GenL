from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genl import Control, Instrument, Layer, calc_kinematic
from genl.paths import FORM_FACTOR_DIR, STRUCTURE_DIR

theta = np.arange(20.0, 90.0, 0.02)
wavelength = 1.54056
q = 4.0 * np.pi / wavelength * np.sin(np.deg2rad(theta))

layer = Layer(
    direction=3,
    n=40,
    filename="Fe_fractional.vasp",
    scale=1.0,
    area_scale=1.0,
)

result = calc_kinematic(
    q,
    wavelength,
    [layer],
    Control(pol=0),
    Instrument(theta=theta),
    poscar_dir=STRUCTURE_DIR,
    form_factor_dir=FORM_FACTOR_DIR,
)

print(f"calculated {len(result.refl)} intensity points")
print(f"min={result.refl.min():.6e} max={result.refl.max():.6e}")
