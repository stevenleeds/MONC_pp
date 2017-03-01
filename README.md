**NOTE**

This code works with a version of MONC that puts the profiles of reference
variables in the model dumps

**LICENSE**

*These programs/scripts are free software: you can redistribute it and/or modify*
*it under the terms of the GNU Lesser General Public License as published by*
*the Free Software Foundation, either version 3 of the License, or*
*(at your option) any later version.*

*This program is distributed in the hope that it will be useful,*
*but WITHOUT ANY WARRANTY; without even the implied warranty of*
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*
*GNU General Public License for more details.*

*You should have received a copy of the GNU Lesser General Public License*
*along with this program.  If not, see <http://www.gnu.org/licenses/>.*

*Copyright (C) 2013-2014 Steven Boeing, ETHZ*

*Copyright (C) 2015 Steven Boeing, University of Leeds*

**CONTACT**

S."lastname with oe" (at) leeds.ac.uk
OR
sjboing (at) "the usual g-mail suffix"

**DESCRIPTION**

This is a collection of scripts for post-processing (after-burning) MONC output

This directory includes the following scripts
* pp_MONC.py: Generic postprocessing routine for a single case. This is the code
  that users can modify easily
* pp_MONC_infrastucture.py: Classes and functions used in pp_MONC.py
* pp_MONC_constants.py: Physical constants from MONC, for calculation of derived variables
* pp_MONC_settings.py: Set what type of output to generate, and at which locations. 
  Also sets part of the file path on different machines
  Uses .cfg configuration files
* fieldslist.py: Lists of variables to analyze
* areaspectra.py: code for calculating spectra (thanks Juerg Schmidli, ETHZ)
* .syscfg files: the settings on different systems (e.g. paths)
* .cfg files: the settings for different cases (e.g. bubbles on different domains)

**REQUIREMENTS**

Python with netcdf4python, numpy and scipy (for embedding C-code with weave)

**OUTPUT FILES**

* clouds.*.tar: tar file with 3d cloud fields (no gz, as cloud fields themselves are already compressed)
* cross_xy.*.nc: cross-sections in the x,y plane
* cross_xz.*.nc: cross-sections in the x,z plane
* cross_yz.*.nc: cross-sections in the y,z plane
* samp_1d.*.nc: sampled profile diagnostics
* spec_x.*.nc: spectra along the x-direction
* spec_y.*.nc: spectra along the y-direction
* stat_1d.*.nc: mean profile diagnostics
* stat_dom.*.nc: domain integrated diagnostics
* stat_int.*.nc: column integrated diagnostics
* stat_xz.*.nc: mean diagnostics in the x,z plane
* stat_yz.*.nc: mean diagnostics in the y,z plane

**ADVANTAGES**

* Easy to add a variable, and obtain a range of outputs on it.
* Correction for staggered grids already implemented.
* Easy to add e.g. vector operations exploiting standard python 
  (see pp_MONC_infrastucture.py for examples, such as the calculation
  of vortivity in the hvort and vvort routines).
* Optional use of netcdf compression
  
**KNOWN ISSUES/FEATURES**

* The code works on 3D NetCDF files that include reference profiles for pressure and
  density, as well as the default checkpoint outputs. The reference checkpoint
  files are not included in checkpoints on the trunk, though a ticket has
  been raised to address the issue. SB has a version of the code that does
  output checkpoitn profiles.
* For very large cases, one could also run out of memory. Suggestions:
  - use node with more memory
  - delete temporary variables using "del", or restructure script
* The code works in serial. In theory, it would be easy to have different nodes/
  cores could work on different parts of the data (memory permitting). It would
  be possible to restructure the code to use e.g. numba as well, but so far
  there has been no need to.

**EXAMPLE USAGE**

python pp_MONC.py bomex small_100
"bomex" is the case directory name
"small_100" is the experiment name

python pp_MONC.py bomex small_100 -c default.cfg -s xcm.syscdf
"bomex" is the case directory name
"small_100" is the experiment name
"default.cfg" is the configuration
"xcm.syscfg" sets the paths on the xcm

The input files are located at e.g.
/nfs/see-fs-01_users/"yourname"/MONCin/bomex/small_100/ (on the SEE systems in Leeds) 

**AVAILABLE VARIABLES**

* See fieldslist.py for a list of fields and their description.

**ADDING A VARIABLE**

* Append the variable to the corresponding list in fieldslist.py (see the top of this file for more details)

* Add the output function and the calculation of the variable to pp_MONC.py. Output functions are

  process_var: process statistics and make cross sections

  stat_var: process only the statistics of the variable

  int_var: this is a height integrated variable

  ref_var: this is a reference variable profile

**RETRIEVING VARIABLES FROM INPUT NETCDF**

* The following functions are used: 

  self.gv: get a variable
  
  self.gq: get a scalar
  
  self.gref: get a reference profile

**CHANGING THE PATH OF INPUT/OUTPUT**

* See pp_MONC_settings.py for default paths on a number of platforms
