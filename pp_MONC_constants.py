# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# Copyright (C) 2013-2015 Steven Boeing, ETHZ
# Copyright (C) 2015 Steven Boeing, University of Leeds

# A PYTHON SCRIPT TO POSTPROCESS ALL MONC OUTPUT IN A DIRECTORY AT ONCE
# BASED ON EARLIER WORK BY THE AUTHOR FOR COSMO OUTPUT

# SEPARATE SCRIPTS SHOULD DO THE PLOTTING OF OUTPUT
# AIMED AT A COMBINATION OF SPEED; READABILITY AND A SOMEWHAT LIMIT MEMORY USAGE
# I.E. ABILITY TO POSTPROCESS LARGE DATA ON A SINGLE NODE

# detect all MONC output files and find the corresponding time
# incrementally add the output to an output file
# includes masked/sampled statistics (like in traditional LES models)
# masks can also be used to average over pre-selected subdomain (e.g. a box)
# but currently such domains are not implemented
# outputs to a zipped file which is saved with f4 precision (to reduce file size)
# and further reduced precision for in-cloud 3d variables

# Requires python with netcdf4python, numpy and scipy (for embedding C-code with weave)

# EXAMPLE
# python pp_MONC.py bomex bomex
# ------------------case directory
# ------------------------experiment name                        

# Moisture vairables
# Negative numbers can be used to ignore species
nqv=0  # q number corresponding to water vapor
nqc=1  # q number corresponding to cloud droplet liquid water
nqi=-1 # q number corresponding to cloud ice
nqs=-1 # q number corresponding to snow
nqg=-1 # q number corresponding to graupel
nqh=-1 # q number corresponding to hail
   
##  MONC CONSTANTS FOR DERIVED VARIABLES

psfr=1.0e5 # reference pressure
rd=287.05 # gas constant, dry air
rvord=1.608 # ratio of gas constants (water vapor/dry air)
cp=1005.0 # heat capacity of air at constant pressure
rlvap=2.501e6 # latent heat of condensation
rlsub=2.834e6 # latent heat of sublimation
grav=9.81 # gravitational acceleration

# Constants used in calculation of saturation pressure (MONC specific)
# qsat=qsa1/(p*exp(qsa2*(t - tk0c)/(T - qsa3)) - qsa4) 
tk0c = 273.15      # Temperature of freezing in Kelvin
qsa1 = 3.8         # Numerator in equation to calculate qsat
qsa2 = -17.2693882 # Constant in qsat equation
qsa3 = 35.86       # Constant in qsat equation
qsa4 = 6.109       # Constant in qsat equation 
qis1 = 3.8         # Numerator in equation to calculate qsat
qis2 = -21.8745584 # Constant in qisat equation
qis3 = 7.66        # Constant in qisat equation
qis4 = 6.109       # Constant in qisat equation
