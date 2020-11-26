import numpy as np
from voroana.voronoi import *
from voroana.utils import cell_to_hmatrix
from ase.data import atomic_numbers
from ase.geometry import find_mic


__all__ = ["cell_to_hmatrix", "get_atomic_radius", "voronoi_analysis"]


def get_atomic_radius(symbols):
    # -----------------------------------------------------------------------------
    #  Atomic radius data taken from:
    #  Wikipage of Atomic radius of elements
    # 
    #  As atomic radius are not uniquely determined, the priority in selecting
    #  the value for each element is as follows:
    #    Metals    : Metallic radius -> empirical -> van der Waals -> Calculated
    #    Non-metal : Covalent -> empirical -> van der Waals -> Calculated
    #    RareGas   : van der Waals -> empirical -> Calculated -> Covalent
    #    NoData    : 1.111
    # ---------------------------------------------------------------------------*/
    atomic_radius = \
    [
        1.111,                             # Null
        0.380, 0.320, 1.520, 1.120, 0.820, # H  - B   1  - 5
        0.770, 0.750, 0.730, 0.710, 1.540, # C  - Ne  6  - 10
        1.860, 1.600, 1.430, 1.110, 1.060, # Na - P   11 - 15
        1.020, 0.990, 1.880, 2.270, 1.970, # S  - Ca  16 - 20
        1.620, 1.470, 1.340, 1.280, 1.270, # Sc - Mn  21 - 25
        1.260, 1.250, 1.240, 1.280, 1.340, # Fe - Zn  26 - 30
        1.350, 1.220, 1.190, 1.160, 1.140, # Ga - Br  31 - 35
        2.020, 2.480, 2.150, 1.800, 1.600, # Kr - Zr  36 - 40
        1.460, 1.390, 1.360, 1.340, 1.340, # Nb - Rh  41 - 45
        1.370, 1.440, 1.510, 1.670, 1.450, # Pd - Sn  46 - 50
        1.450, 1.400, 1.330, 2.160, 2.650, # Sb - Cs  51 - 55
        2.220, 1.870, 1.818, 1.824, 1.814, # Ba - Nd  56 - 60
        1.834, 1.804, 1.804, 1.804, 1.773, # Pm - Tb  61 - 65
        1.781, 1.762, 1.761, 1.759, 1.760, # Dy - Yb  66 - 70
        1.738, 1.590, 1.460, 1.390, 1.370, # Lu - Re  71 - 75
        1.350, 1.355, 1.385, 1.440, 1.510, # Os - Hg  76 - 80
        1.700, 1.800, 0.160, 1.900, 2.020, # Tl - At  81 - 85
        2.200, 3.480, 2.150, 1.950, 1.790, # Rn - Th  86 - 90
        1.630, 1.560, 1.550, 1.590, 1.730, # Pa - Am  91 - 95
        1.740, 1.700, 1.860, 1.860, 1.111, # Cm - Fm  95 - 100
        1.111, 1.111, 1.111, 1.310, 1.260, # Md - Db  101- 105
        1.210, 1.190, 1.180, 1.130, 1.120, # Sg - Ds  106- 110
        1.180, 1.300                       # Rg - Cn  111- 112
    ]
    
    return {i: atomic_radius[atomic_numbers[i]] for i in set(symbols)}


def guess_grid(n, lx, ly, lz):
    l = (n/(5.6*lx*ly*lz))**(1/3)
    nx = int(lx*l + 1)
    ny = int(ly*l + 1)
    nz = int(lz*l + 1)

    return nx, ny, nz


def refine_cell(index, cell, atoms, cutoff, relative):
    lx, _, ly, _, _, lz = a.cell[np.tril_indices(3)].tolist()
    fa = cell.face_areas()
    fs = 1.0
    if relative:
        fs = 1/sum(fa)
    fc = [i*fs > ratio for i in fa]
    if all(fc):
        return cell
    else:
        cnew = Cell()
        cnew.init(-lx, lx, -ly, ly, -lz, lz)
        neighbors = cell.neighbors()
        normals = np.array(cell.normals())
        _, d = find_mic(atoms.positions[neighbors, :] - atoms.positions[index, :], atoms.cell)
        vectors = normals[:,:]*d[:,np.newaxis]
        for i in range(len(neighbors)):
            if neighbors[i] < 0:
                continue
            if fc[i]:
                assert cnew.nplane(vectors[i, 0], vectors[i, 1], vectors[i, 2], neighbors[i])
        return cnew


def voronoi_analysis(atoms, outputs="ln", max_face_orders=6, radius=None, cutoff=0.02, relative=True):
    if not atoms.pbc.all():
        raise RuntimeError("Only full pbc conditions supported so far!")
    a = cell_to_hmatrix(atoms)
    bx, bxy, by, bxz, byz, bz = a.cell[np.tril_indices(3)].tolist()
    nx, ny, nz = guess_grid(len(a), bx, by, bz)
    if np.allclose(a.cell[np.tril_indices(3, k=-1)], 0.0):
        if radius == None:
            con = Container(0, bx, 0, by, 0, bz, nx, ny, nz, True, True, True, 8)
            for i in a:
                con.put(i.index, i.x, i.y, i.z)
        else:
            con = ContainerPoly(0, bx, 0, by, 0, bz, nx, ny, nz, True, True, True, 8)
            for i in a:
                con.put(i.index, i.x, i.y, i.z, radius[i.symbol])
    else:
        if radius == None:
            con = ContainerPeriodic(bx, bxy, by, bxz, byz, bz, nx, ny, nz, 8)
            for i in a:
                con.put(i.index, i.x, i.y, i.z)
        else:
            con = ContainerPeriodicPoly(bx, bxy, by, bxz, byz, bz, nx, ny, nz, 8)
            for i in a:
                con.put(i.index, i.x, i.y, i.z, radius[i.symbol])
    cells = con.get_cells()
    assert len(cells) == len(a)
    for i in range(len(cells)):
        cells[i] = refine_cell(cells[i].id, cells[i], a, cutoff, relative)

    names = {
        "i": "indices",
        "l": "labels",
        "n": "neighbors",
        "N": "neighbor_vectors",
        "s": "face_areas",
        "v": "volumes"
    }

    results = {names[i]: [] for i in outputs}
    for cell in cells:
        t1 = cell.face_orders()
        t2 = [t1.count(i) for i in range(3, max_face_orders + 1)]
        for i in outputs:
            if i == "i":
                results[names[i]].append(t2)
            if i == "l":
                results[names[i]].append(",".join([str(i) for i in t2]))
            if i == "n":
                results[names[i]].append(cell.neighbors())
            if i == "N":
                neighbors = cell.neighbors()
                normals = np.array(cell.normals())
                _, d = find_mic(a.positions[neighbors, :] - a.positions[cell.id, :], a.cell)
                vectors = normals[:,:]*d[:,np.newaxis]
                results[names[i]].append(vectors)
            if i == "s":
                results[names[i]].append(cell.face_areas())
            if i == "v":
                results[names[i]].append(cell.volume())

    return results

