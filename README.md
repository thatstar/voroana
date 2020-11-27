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
   

Absolute or relative face area based cutoff can be applied to remove tiny cell faces during the calculation and radical Voronoi tessellations can be performed on requirement.



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

>`ase >= 3.20.1`
>
>`python >= 3.7`
>
>(optional) `matplotlib > 3.0`

## 3. Running analysis

```python
from ase.build import bulk
from voroana import get_atomic_radius, voronoi_analysis

a = bulk("Pd")
v = voronoi_analysis(a,
                     outputs="ilnsv",
                     radius=get_atomic_radius(a.get_chemical_symbols()),
                     max_face_order=6,
                     cutoff=0.0,
                     relative=True)
for i in range(len(a)):
    print(a[i].symbol, v['labels'][i])
```

For additional information,  open an `IPython` console and type `:help voronoi_analysis` or consult the source codes of the package.



# Known limitations

Currently, only Atoms object with full PBC condition (`atoms.pbc = True`) is supported.



# License

The code is released under the GPLv3 license while the voro++ code is released under the BSD license.



