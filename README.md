# grid_continuity_SMB
Collection of MATLAB codes and scripts to estimate SMB based on remote sensing observations of velocity and surface elevation changes, with ice thickness estimates. 

## Prerequisites

You need:
 - activated version of MATLAB (R2017a or more recent)
 - Mapping Toolbox to ingest Geotiff information

Other dependencies are included in /code

## Installing

Simply clone the repository, then adjust the master_FluxCalcsSimple.m script to get going. Example data are included for Matanuska and Rongbuk Glaciers. Ensure that input velocity data (dx, dy) are relative to their coordinate system.

# Authors

* **Evan Miles** - *Flux-based SMB calculations* - [miles916](https://github.com/miles916), [Miles_of_Ice](https://twitter.com/Miles_of_Ice)
* **Amaury Dehecq** - *Conceptual development* - [adehecq](https://github.com/adehecq)

# License

Currently under private development.

# Acknowledgments
*Thanks to [ColorBrewer](http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3) for their development of improved color ramps; here implemented with a [FileExchange script](https://ch.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab).
