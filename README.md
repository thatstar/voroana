Voroana, a python package for calculating the Voronoi statistics of atom configurations.

by Su Rui, surui@hdu.edu.cn



[TOC]

#  About

The **voroana** library consists as a cython interface to the **voro++** library and high-level python functions for calculating the Voronoi based statistics for atomistic simulations. Currently, it supports to calculate the following properties:

1. Voronoi indices up to user defined maximum face order.
2. Voronoi (atomic) volumes.
3. Neighbor list based on the Voronoi tessellation.
4. Neighbor vectors.
5. Face areas of the Voronoi cells.
6. Vectors from atoms to their cell centroids.

Absolute or relative face area based cutoff can be applied to remove tiny cell faces during the calculation and radical Voronoi tessellations can be performed on requirement.

Additionaly, a set of analysis functions based on Voronoi analysis are provided:

1. Bond orient order ($q_{l}$, $Q_{l}$, $W_{l}$, $\hat{W}_{l}$) for arbitrary angular moment order $l$  with the neighbor information from Voronoi analysis. Face area weighted version can be specified if required. 
2. self intermediate scattering function for all or selected atoms with a specified voronoi index.




# Usage

## 1. Installation

The package can be installed in the Linux environment by following procedure:

```bash
tar xvf voroana-xxx.tar.gz
cd voroana-xxx
python setup.py build_ext --inplace
```

Then add a line:`export PYTHONPATH=/path/to/voroana-xxx:$PYTHONPATH` in your `.bashrc` file and `source ~/.bashrc`.

## 2. Dependence

>`python >= 3.7`
>
>`ase >= 3.20.1`
>
>`sympy >= 1.2`
>
>`scipy >= 0.9`
>
>(optional) `matplotlib > 3.0`
>
>(optional) `seaborn > 0.10`

## 3. Running analysis

```python
from ase.build import bulk
from voroana import get_atomic_radius, voronoi_analysis
from voroana.analysis import BondOrientOrder

a = bulk("Pd", cubic=True)
v = voronoi_analysis(a,
                     outputs="ilnsv",
                     radius=get_atomic_radius(a.get_chemical_symbols()),
                     max_face_order=6,
                     cutoff=0.0,
                     relative=True)
for i in range(len(a)):
    print(a[i].symbol, v['labels'][i])

boo = BondOrientOrder(a, voronoi_analysis(a, outputs="snN"), weighted=True)
print(boo.compute_q_l(4, coarse_grained=False)) # q4
print(boo.compute_q_l(4, coarse_grained=True)) # Q4
print(boo.compute_w_l(6, coarse_grained=True)) # W6
```

For additional information,  open an `IPython` console and type `:help voronoi_analysis` or consult the source codes of the package.



# Known limitations

1. Currently, only Atoms object with full PBC condition (`atoms.pbc = True`) is supported.
2. The cell must be large enough so that the images of the center atom do not present in its neighbors.



# License

The code is released under the GPLv3 license while the voro++ code is released under the BSD license.



