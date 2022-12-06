# xsphere - A package to display NWP fields of unstructured grid with xarray.

The code in this repository provides a solution to display with xarray NWP fields of unstructured grid. 
It builds mainly upon xarray, cartopy and pygsp.

ATTENTION: The code is subject to changes in the coming months.

The folder `tutorials` (will) provide jupyter notebooks describing various features of xsphere.

The folder `docs` (will) contains slides and notebooks explaining the xsphere framework.

## Installation

For a local installation, follow the below instructions.

1. Clone this repository.
   ```sh
   git clone https://github.com/ghiggi/xsphere.git
   cd xsphere
   ```

2. Install manually the following dependencies:
   ```sh
   conda create --name xsphere-dev python=3.8
   conda install xarray h5netcdf netcdf4 zarr 
   conda install cartopy matplotlib-base
   conda install scipy numpy
   pip install git+https://github.com/epfl-lts2/pygsp@sphere-graphs
   ```
   
2. Alternatively install the dependencies using one of the appropriate below 
   environment.yml files:
   ```sh
   conda env create -f TODO.yml
   ```

## Tutorials


## Contributors

* [Gionata Ghiggi](https://people.epfl.ch/gionata.ghiggi)
 

## License

The content of this repository is released under the terms of the [MIT license](LICENSE.txt).
