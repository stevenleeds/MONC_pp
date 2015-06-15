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
* fieldslist.py: Lists of variables to analyze
* areaspectra.py: code for calculating spectra (thanks Juerg Schmidli, ETHZ)

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

**KNOWN ISSUES/FEATURES**

* For very large cases, one could also run out of memory (use node with more memory,
  delete temporary variables using "del", or restructure script) or processing time

**EXAMPLE USAGE**

python pp_MONC.py bomex bomex
The first "bomex" is the case directory name
The second "bomex" is the experiment name

With the input files in e.g.
/nfs/see-fs-01_users/"yourname"/MONCin/bomex (at Leeds)

**ADDING A VARIABLE**

* Append the variable to the corresponding list in fieldslist.py
* Add the output function and the calculation of the variable to pp_MONC.py. Output funtions are
**  process_var: process statistics and make cross sections
**  stat_var: process only the statistics of the variable
**  int_var: this is a height integrated variable
**  ref_var: this is a reference variable profile

**CHANGING THE PATH OF INPUT/OUTPUT**

* See pp_MONC_settings.py for default paths on a number of default paths
