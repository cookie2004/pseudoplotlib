# PseudoPlotLib
# Forked from @sokrypton

<img src='./docs/idp_traj-ezgif.com-crop.gif' width='600' style="horizontal-align:middle">

This repository provides tools to generate pseudo protein structures.  

## Installation

You can install the package directly from GitHub:

```bash
pip install git+https://github.com/cookie2004/pseudoplotlib.git
```

Or clone the repository and install locally:

```bash
git clone https://github.com/cookie2004/pseudoplotlib.git
cd pseudoplotlib
pip install -e .
```

## Requirements

- Python >= 3.6
- NumPy >= 1.18
- Matplotlib >= 3.2

## Features

PseudoPlotLib is a specialized plotting library designed for visualizing pseudo protein representations. Built on top of matplotlib, it offers intuitive functions for creating detailed visualizations of protein structures and related data.

## Usage

```python
import pseudoplotlib.pyplot as pplt
import pseudoplotlib.utils as pplutils
import matplotlib.pylab as plt

# Basic example
# Load and convert PDB file to XYZ coordinates
monomer_filename = "tests/monomer/1clm.pdb"
xyz_coords = pplutils.pdb_to_xyz(monomer_filename)
# Plot the XYZ coordinates
fig, ax = plt.subplots()
#cmap can be your favorite matplotlib cmap!
pplt.pseudo(xyz_coords, ax=ax, cmap='viridis')
plt.show()

# Creating animation
# Load and convert PDB file to XYZ coordinates
idp_traj_filename = "tests/traj/idp_traj.pdb"
xyz_coords = pplt.pdb_traj_to_xyz(idp_traj_filename)
# Plot the XYZ coordinates
animation = pplt.make_animation(xyz_coords, ax=ax)
animation.save('tests/traj/idp_traj.gif')

```

## Documentation

More detailed documentation and examples will be added as the project develops.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

- Author: KittenMittens
- Email: cookie2004@gmail.com
- GitHub: [cookie2004](https://github.com/cookie2004)
