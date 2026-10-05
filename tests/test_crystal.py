from pathlib import Path
import tempfile

import numpy as np

from genl.crystal import project_structure, write_genl_projection
from genl.poscar import PoscarStructure, read_poscar
from genl.paths import STRUCTURE_DIR


def test_cartesian_selective_poscar_is_normalized_to_fractional():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "POSCAR"
        path.write_text(
            "X\n2.0\n1 0 0\n0 1 0\n0 0 1\nSi\n1\n"
            "Selective dynamics\nCartesian\n"
            "0.5 1.0 1.5 T F T\n",
            encoding="utf-8",
        )
        structure = read_poscar(path)
        np.testing.assert_allclose(structure.lattice, np.eye(3) * 2.0)
        np.testing.assert_allclose(structure.positions[0], [0.5, 1.0, 1.5])


def test_si_111_projection_matches_known_bilayer_and_existing_oriented_cell():
    conventional = read_poscar("Si_fractional.vasp", STRUCTURE_DIR)
    projected = project_structure(conventional, (1, 1, 1))

    d111 = 5.431 / np.sqrt(3.0)
    assert np.isclose(projected.period, d111)

    levels, counts = np.unique(
        np.round(projected.z / projected.period, 10),
        return_counts=True,
    )
    np.testing.assert_allclose(levels, [0.0, 0.75])
    np.testing.assert_array_equal(counts, [4, 4])

    oriented = read_poscar("si_111_fractional.vasp", STRUCTURE_DIR)
    oriented_z = oriented.positions @ oriented.a3
    oriented_area = np.linalg.norm(np.cross(oriented.a1, oriented.a2))

    # Compare normalized specular structure factors at several (111) orders.
    for order in range(1, 7):
        q = order * 2.0 * np.pi / projected.period
        projected_f = np.sum(np.exp(1j * q * projected.z)) / projected.area
        oriented_f = np.sum(np.exp(1j * q * oriented_z)) / oriented_area
        np.testing.assert_allclose(projected_f, oriented_f, rtol=1e-7, atol=1e-9)


def test_projection_writer_round_trips_genl_geometry():
    structure = PoscarStructure(
        ("X",),
        np.array([1]),
        np.array([[0.0, 0.0, 0.0]]),
        *np.eye(3) * 3.0,
    )
    projected = project_structure(structure, (1, 1, 1))

    with tempfile.TemporaryDirectory() as directory:
        path = write_genl_projection(
            projected,
            Path(directory) / "projected.vasp",
        )
        reread = read_poscar(path)
        assert np.isclose(np.linalg.norm(reread.a3), projected.period)
        assert np.isclose(
            np.linalg.norm(np.cross(reread.a1, reread.a2)),
            projected.area,
        )
