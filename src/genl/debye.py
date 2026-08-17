from __future__ import annotations

import numpy as np

_DEBYE_WALLER_B = {
    3: 4.1,
    11: 7.9,
    13: 0.86,
    14: 0.45,
    19: 12.0,
    23: 0.55,
    24: 0.26,
    26: 0.35,
    28: 0.37,
    29: 0.57,
    32: 0.57,
    41: 0.49,
    42: 0.25,
    46: 0.45,
    47: 0.79,
    73: 0.32,
    74: 0.18,
    78: 0.32,
    79: 0.57,
    82: 2.42,
}


def debye_waller_prefactor(z: int | np.ndarray) -> float:
    """Return Sears room-temperature B factor used by MATLAB `DB_prefactor`."""

    z_values = np.atleast_1d(np.asarray(z, dtype=int))
    if len(z_values) != 1:
        return 0.0
    return _DEBYE_WALLER_B.get(int(z_values[0]), 0.0)
