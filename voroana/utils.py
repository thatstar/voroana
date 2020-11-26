import numpy as np
from ase import Atoms
from ase.geometry import wrap_positions


def cell_to_hmatrix(atoms):
    """
    Transform an ASE atoms object into a regularized one comforming the LAMMPS
    specifications.

    Arguments:
        atoms: ASE atoms object

    Returns:
        Transformed atoms.
    """
    H, pos = cell_to_hmatrix_impl(atoms.cell, atoms.positions, atoms.pbc)

    atoms_new = Atoms(atoms.get_chemical_symbols(), cell=H, positions=pos)

    return atoms_new


def cell_to_hmatrix_impl(cell, positions, pbc):
    """
    Rotate ASE atoms cell to the h-matrix format used by LAMMPS.

    Arguments:
        cell: 3x3 array-like, the 3x3 cell matrix.
        positions: array-like, positions of atoms.
        pbc: tuple of bools, whether the cell is periodic in the direction?

    Returns:
        cell and positions under the new h-matrix coordination system.
    """
    a, b, c = cell
    an, bn, cn = [np.linalg.norm(v) for v in cell]

    alpha = np.arccos(np.dot(b, c)/(bn*cn))
    beta  = np.arccos(np.dot(a, c)/(an*cn))
    gamma = np.arccos(np.dot(a, b)/(an*bn))

    xhi = an
    xyp = np.cos(gamma)*bn
    yhi = np.sin(gamma)*bn
    xzp = np.cos(beta)*cn
    yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
    zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    Hpre = np.array([[xhi, 0.0, 0.0], [xyp, yhi, 0.0], [xzp, yzp, zhi]])
    R = np.dot(np.linalg.inv(cell), Hpre)

    def fold(vec, pvec, i):
        p = pvec[i]
        x = vec[i] + 0.5*p
        n = (np.mod(x, p) - x) / p
        return vec + n*pvec

    Hpre[1, :] = fold(Hpre[1, :], Hpre[0, :], 0)
    Hpre[2, :] = fold(Hpre[2, :], Hpre[1, :], 1)
    Hpre[2, :] = fold(Hpre[2, :], Hpre[0, :], 0)

    H = Hpre

    rot_positions = np.dot(positions, R)

    rot_positions = wrap_positions(rot_positions, H, pbc)

    return H, rot_positions


if __name__ == "__main__":
    from ase.visualize import view
    from ase.build import bulk

    au = bulk("Au") * (3, 3, 3)

    au_new = cell_to_hmatrix(au)

    view(au_new)

