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

import os
import socket
import getpass

## SETTINGS

loadmpl=False # Load matplotlib
lzlib=True # Compress output using zlib?
lcross_xy=True # Produce xy cross-sections? Mais oui
lcross_xz=True # Produce xz cross-sections? Mais oui
lcross_yz=True # Produce yz cross-sections? Mais oui
ldiag_xz=True # Produce xy diagnostics? Mais oui
ldiag_yz=True # Produce yz diagnostics? Mais oui
ldiag_int=True # Produce vertically integrated diagnostics? Mais oui
lsamp=True # Produce sampled diagnostics? Mais oui
lspec=True # Produce spectral diagnostics? Mais oui
lclouds=True # Produce tar-ball with 3D cloud and precipitation fields for storage? Mais oui

# where to take cross sections (currently just grid numbers, i.e., no interpolation)
# has to be an array
xsel=[0]
ysel=[0]
zsel=[0,1,10,20,30]

# Between which level to produce spectra (can be multiple sets, formatted as [lowerlevel,upperlevel+1], as python array indexing works)
spectralevels=[
[10,20],
[20,30],
]

# Option to crop the boundaries in the horizontal plane
nboundlines=0

## FILE PATHS

myusername=getpass.getuser()
hostname=socket.gethostname()
homedir = os.environ['HOME']

# fullbase: standard location of input directories
# scratchbase: location of scratch space where initial postprocessing is done
# localbase: location where final files are stored

if 'see' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
    scratchbase='/scratch/MONCout/'
    projectbase='/nfs/see-fs-01_users/'+myusername+'/MONCout/'
    fullbase='/nfs/see-fs-01_users/'+myusername+'/MONCin/'
elif 'arc2' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
    scratchbase='/nobackup/'+myusername+'/MONCout/'
    projectbase='/home/ufaserv1_h/'+myusername+'/MONCout/'
    fullbase='/nobackup/'+myusername+'/'
elif 'monsoon-metoffice' in hostname: ## Met Office monsoon
    scratchbase='/scratch/'+myusername+'/MONCout/'
    projectbase='/projects/udmonc/'+myusername+'/MONCout/'
    fullbase='/scratch/'+myusername+'/'
else: ## e.g. laptop
    fullbase=homedir+'/MONCin/'
    scratchbase=homedir+'/scratch/MONCout/'
    projectbase=homedir+'/proj/MONCout/' 
