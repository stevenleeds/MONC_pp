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
# and further reduced precision for in-cloud variables

# Requires python with netcdf4python, numpy and scipy (for embedding C-code with weave)

# EXAMPLE
# python pp_MONC.py bomex bomex
# ------------------case directory
# ------------------------experiment name                        

# TODO
# ADD FURTHER VARIABLES (NEEDS REFERENCE PROFILES)
# ADD SAMPLING BASED ON BUOYANCY
# SAMPLING OF STAGGERED VARIABLES...CURRENTLY USING DIFFERENT GRIDS FOR MASKS
#   INTERPOLATE SATURATION DEFICIT AND DETERMINE BUOYANCY USING INTERPOLATED CONSERVED VARIABLES?
# FIGURE OUT Q NUMBER FOR ICE (AND FOR OTHER SPECIES)
# DISCUSS MOMENTUM, ENERGY, VORTICITY, HELICITY AND "PRESSURE LAPLACIAN" BUDGET TOOLS
# ADD AN OPTION PARSER?
# TAKE STAGGERING/WEIGHTING INTO ACCOUNT TO PRODUCE SPECTRA AT THE SAME LEVELS FOR W AND OTHER VARIABLES?

from numpy import *
from netCDF4 import Dataset
from scipy import weave # to embed code in C
from scipy.weave import converters
from fieldslist import * # a separate file contains the actual variable list
from contextlib import closing
import sys # system functions
import glob # a libary for regular expressions in file names
import shutil # for copying the files to destination
import areaspectra # a separate library for computing spectra
import tarfile # for compressing the cloud field data
import os
import errno
import time
import getpass
import numpy
import socket
   
nboundlines=0 # Option to crop the boundaries in the horizontal plane

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

nqc=1 # q number corresponding to cloud droplet liquid water
nqi=-1 # q number corresponding to cloud ice

myusername=getpass.getuser()
hostname=socket.gethostname()
homedir = os.environ['HOME']

## STORAGE LOCATIONS
# fullbase: standard location of input directories
# scratchbase: location of scratch space where initial postprocessing is done
# localbase: location where final files are stored

if 'see' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
    scratchbase='/scratch/MONCout/'
    projectbase='/nfs/see-fs-01_users/'+myusername+'/MONCout/'
    fullbase='/nfs/see-fs-01_users/'+myusername+'/MONCin/'
if 'arc2' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
    scratchbase='/nobackup/'+myusername+'/'
    projectbase='/home/ufaserv1_h/'+myusername+'/MONCout/'
    fullbase='/nobackup/'+myusername+'/'
else: ## E.G. LAPTOP
    fullbase=homedir+'/MONCin/'
    scratchbase=homedir+'/scratch/MONCout/'
    projectbase=homedir+'/proj/MONCout/'
    
# Some MONC constants for calculations of derived variables
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

start=time.clock()
numpy.seterr(invalid='ignore') # don't wine about nans

# currently not used, but may be useful when analysing problems
if loadmpl:
    import matplotlib
    matplotlib.use('agg') # first define agg output, then continue to load rest of matplotlib and pyplot 
    from matplotlib import *
    from matplotlib.pyplot import *
    
# forced makedir
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# make a filelist of the 3D output
def make_filelist():
    global filelist
    types = {'3ddump':'*.nc'} # file type list, currently contains only 3D output to be postprocessed
    filelist={}
    for i in types:
        filelist[i]=glob.glob(fulldir+types[i])
        filelist[i].sort() 
        
# create a tarfile (for cloud fields)
def make_tarfile(output_filename, infiles):
    with closing(tarfile.open(output_filename,'w')) as tar:
        for infile in infiles:
            tar.add(infile)
                          
# class for masks/conditional sampling
class mask:
    def __init__(self):
        pass
    def setfield(self,field):
        self.field=field
    def setfieldxe(self,field):
        self.fieldxe=field
    def setfieldye(self,field):
        self.fieldye=field
    def setfieldze(self,field):
        self.fieldze=field
                                             
# prepare spectral data into 2d arrays
# containing stretches along all points 
# that match sampling criteria
# this is for y-direction geometry
# uses C-code (weave) for speed
def prepare_spectra_y(dataout,data,imax,jmax,kmin,kmax):
    code = """
    int i, j, k, nsegments;
    nsegments=0;
    for (i=0; i<imax; i++) {
        for (k=kmin; k<kmax; k++) {
            for (j=0; j<jmax; j++) {
                dataout(nsegments,j)=data(i,j,k);
            }
            nsegments=nsegments+1;
        }     
    }
    """
    weave.inline(code,['dataout','data','imax','jmax','kmin','kmax'],type_converters=converters.blitz,compiler='gcc')

# The spectra with along x-direction
def prepare_spectra_x(dataout,data,imax,jmax,kmin,kmax):
    code = """
    int i, j, k, nsegments;
    nsegments=0;
    for (j=0; j<jmax; j++) {
        for (k=kmin; k<kmax; k++) {
            for (i=0; i<imax; i++) {
                dataout(nsegments,i)=data(i,j,k);
            }
            nsegments=nsegments+1;
        }     
    }
    """
    weave.inline(code,['dataout','data','imax','jmax','kmin','kmax'],type_converters=converters.blitz,compiler='gcc')
             
# mean value across 2d slab
def mean_2d(inputfield):
    meanfield=mean(mean(inputfield,axis=1, dtype=numpy.float64),axis=0, dtype=numpy.float64)
    return meanfield
    
# mean value across 1d line      
def mean_xz(inputfield):
    meanfield=mean(inputfield,axis=1, dtype=numpy.float64)
    return meanfield
    
def mean_yz(inputfield):
    meanfield=mean(inputfield,axis=0, dtype=numpy.float64)
    return meanfield

# extraction used for cross-sections
def extract_xz(inputfield,indices):
    extrfield=inputfield[:,indices,:]
    return extrfield
    
def extract_yz(inputfield,indices):    
    extrfield=inputfield[indices,:,:]
    return extrfield

def extract_xy(inputfield,indices):    
    extrfield=inputfield[:,:,indices]
    return extrfield
       
# deviations across 2d slab   
def deviation_2d(inputfield):
    meanfield=mean_2d(inputfield)
    delta=inputfield-meanfield[None,None,:]
    return delta
                 
# deviations with respect to xz plane     
def deviation_xz(inputfield):
    meanfield=mean_xz(inputfield)
    delta=inputfield-meanfield[:,None,:]
    return delta
    
# deviations with respect to yz plane        
def deviation_yz(inputfield):
    meanfield=mean_yz(inputfield)
    delta=inputfield-meanfield[None,:,:]
    return delta

# command to get a single variable from a file
def var_from_file(dataset,key):
    if(nboundlines>0):
        try:
            return dataset.variables[(key)][nboundlines:-nboundlines,nboundlines:-nboundlines,:]
        except:
            return(nan)
    else:
        try:
            return dataset.variables[(key)][:,:,:]
        except:
            return(nan)
            
# a class to include some methods needed for both the help-variable class (which store temp variables) and the netcdf class
class get_variable_class():
    # gv being "get variable"
    def gv(self,key):
        if(nboundlines>0):
            if key in self.varkeys:
                return self.data.variables[(key)][nboundlines:-nboundlines,nboundlines:-nboundlines,:]
            else:
                return(self.data.variables[('p')][nboundlines:-nboundlines,nboundlines:-nboundlines,:]*nan)
        else:
            if key in self.varkeys:
                return self.data.variables[(key)][:,:,:]
            else:
                return(self.data.variables[('p')][:,:,:]*nan)          
    # gref being "get reference variable"
    def gref(self,key):
        if key in self.varkeys:
            return self.data.variables[(key)][:]
        else:
            return(self.data.variables[('p')][:]*nan) 
    # gq "get moisture variable"
    def gq(self,key,index):
        if(nboundlines>0):
            if key in self.varkeys:
                if index<len(self.data.variables[(key)]):
                    return self.data.variables[(key)][index,nboundlines:-nboundlines,nboundlines:-nboundlines,:]
                else:
                    return zeros(shape(self.data.variables[(key)][0]))
            else:
                return(self.data.variables[('p')][nboundlines:-nboundlines,nboundlines:-nboundlines,:]*nan)
        else:
            if key in self.varkeys:
                return self.data.variables[(key)][index,:,:,:]
            else:
                return(self.data.variables[('p')][:,:,:]*nan)    
    # gdim being "get dimension"
    def gdim(self,key):
        try:
            return self.data.variables[(key)][:]
        except:
            try:
                return xrange(len(self.data.dimensions[(key)]))
            except:
                return([nan]) 
                
# a class to store derived variables from output which are needed relatively often       
class nchelper(object,get_variable_class):
    def __init__(self):
        self.tstep=0
        self.data=[]
        self.varkeys=[]        
    def update(self,data):
        self.data=data
        self.varkeys=self.data.variables.keys()
        self.tstep=self.tstep+1
        w=self.gv('w')
        wplus=dstack((w[:,:,1:],0.0*w[:,:,-1]))
        self.wzc=0.5*(w+wplus)
        self.wmin=nanmin(w,axis=2)
        self.wmax=nanmax(w,axis=2)
        qci=self.gq('q',nqc)+self.gq('q',nqi) # try to include ice
        self.cld=(qci>1.0e-6)
        self.cloudycolumn=1.0*(sum(self.cld,axis=2)>0)
        zmin=self.gdim('z')[:-1]
        zplus=self.gdim('z')[1:]
        zhalf=0.5*(zmin+zplus)
        bottom=-zhalf[0]
        self.zc=hstack(([bottom],zhalf))
                              
class ncobject(object,get_variable_class):
    # class for writing to netcdf
    def __init__(self,outfile,description,logical=True):
        self.active=logical
        if self.active:
            self.data=[]
            self.outvars={}
            self.ncoutname=scratchdir+outfile
            try:
                os.remove(self.ncoutname)
            except:
                pass
            self.outfile=Dataset(self.ncoutname,'w',format='NETCDF4',zlib=lzlib)
            self.outfile.description=description+' for '+case
            self.outfile.history = 'Created ' + time.ctime(time.time())
            self.outfile.source = 'Created by user ' + myusername
            self.outfile.createDimension('time',0)
            timevar=self.outfile.createVariable('time', 'f8', ('time',),zlib=lzlib)
            setattr(timevar,'longname','time [s]')
            self.myvars=[]
            self.outfile.close()
            self.tstep=0     
    def make_dims(self):
        pass
    def make_var(self):
        pass
    # LAYOUT OF POSTPROCESSING A TIME STEP
    def opener(self,data):
        if self.active:
            self.data=data      
            self.outfile=Dataset(self.ncoutname,'a',format='NETCDF4',zlib=lzlib)
            if(self.tstep==0):
                self.set_dims()
                self.t=0
            else:
                self.t=len(self.outfile.variables['time'])
            self.outfile.variables['time'][self.t]=self.data.variables['time'][0]
    def closer(self):
        if self.active:
            self.tstep=self.tstep+1
            self.outfile.close()
    # functions to initialize dimensions in output
    def set_dims(self):
        pass
    def init_dim(self,dimname,longdimname,dimvalues,sel=None):
        if sel==None:
            self.outfile.createDimension(dimname,len(dimvalues))
            var=self.outfile.createVariable(dimname, 'f8', (dimname,),zlib=lzlib)
            var[:]=dimvalues
        else:
            self.outfile.createDimension(dimname,len(sel))
            var=self.outfile.createVariable(dimname, 'f8', (dimname,),zlib=lzlib)
            var[:]=dimvalues[sel]      
        setattr(var,'longname',longdimname)
    def init_dimxc(self,sel=None):
        xe=self.gdim('x')
        xmin=hstack(([0],self.gdim('x')[:-1]))
        self.xc=0.5*(xmin+xe)
        self.init_dim('xc','x [m]',self.xc,sel)        
    def init_dimyc(self,sel=None):
        ye=self.gdim('y')
        ymin=hstack(([0],self.gdim('y')[:-1]))
        self.yc=0.5*(ymin+ye)
        self.init_dim('yc','y [m]',self.yc,sel)
    def init_dimzc(self,sel=None):                      
        zmin=self.gdim('z')[:-1]
        zplus=self.gdim('z')[1:]
        zhalf=0.5*(zmin+zplus)
        bottom=-zhalf[0]
        self.zc=hstack(([bottom],zhalf))
        self.init_dim('zc','height [m]',self.zc,sel)    
    def init_dimxe(self,sel=None):
        self.xe=self.gdim('x') 
        self.init_dim('xe','x [m] (staggered)',self.xe,sel)
    def init_dimye(self,sel=None):
        self.ye=self.gdim('y')
        self.init_dim('ye','y [m] (staggered)',self.ye,sel)
    def init_dimze(self,sel=None):
        self.ze=self.gdim('z')
        self.init_dim('ze','height [m] (bottom staggered)',self.ze,sel)
    def init_varwithdims(self,var,dimsin,mask=None):
        # initialise a variable with correct dimensions on a staggered grid
        dimsarr=['time']
        if 'z' in dimsin:
           if 'ze' in allfields[var]:
              dimsarr+=['ze']
           else:
              dimsarr+=['zc'] 
        if 'y' in dimsin:
           if 'ye' in allfields[var]:
              dimsarr+=['ye']
           else:
              dimsarr+=['yc']        
        if 'x' in dimsin:
           if 'xe' in allfields[var]:
              dimsarr+=['xe']
           else:
              dimsarr+=['xc']
        dims=tuple(dimsarr)
        units=allfields[var][1]
        if mask==None:
            longname=allfields[var][0]
            self.init_var(var,longname,units,dims)
        else:
            longname=allfields[var][0]+' ('+mask+')'
            self.init_var(var+mask,longname,units,dims)         
    def init_var(self,var,longname,units,dims):
        so=self.outfile.createVariable(var, 'f4', dims,zlib=lzlib)
        so.missing_value = nan
        so.long_name=longname
        so.units=units
    # most low level way to write a field
    def put_var(self,var,field,mask=None):
        if mask==None:
            self.outfile.variables[var][self.t]=transpose(field)
        else:
            self.outfile.variables[var+mask][self.t]=transpose(field)
    def make_var(self,var):
        pass # defined in subclasses       
    def put_make_var(self,var,field):
        if self.active:
            if(self.tstep==0):
                self.make_var(var)
            self.put_var(var,field)
    def put_make_sampvar(self,var,field,mask):
        if self.active:
            if(self.tstep==0):
                self.make_var(var,mask)
            temp=1.0*field
            whereinf=isinf(field)
            temp[whereinf] = nan
            self.put_var(var,field,mask)             
    def init_vdims(self,sel=None):
        self.init_dimzc(sel)
        self.init_dimze(sel)
    def init_xdims(self,sel=None):
        self.init_dimxc(sel)
        self.init_dimxe(sel)    
    def init_ydims(self,sel=None):
        self.init_dimyc(sel)
        self.init_dimye(sel)       
                
# general class for 1d statistics based on 3d output
class statgroup_1d(ncobject):    
    def __init__(self,outfile,description):
        super(statgroup_1d,self).__init__(outfile,description)    
    def set_dims(self):
        self.init_vdims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,'z',mask)

# general class for height integrated (min/max/liquid water path etc) output
class statgroup_int(ncobject):    
    def __init__(self,outfile,description):
        super(statgroup_int,self).__init__(outfile,description,ldiag_int)    
    def set_dims(self):
        self.init_hdims()
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims()
    def make_var(self,var):
        self.init_varwithdims(var,['y','x'])

# general class for domain integrated output
class statgroup_dom(ncobject):    
    def __init__(self,outfile,description):
        super(statgroup_dom,self).__init__(outfile,description,ldiag_int)    
    def set_dims(self):
        pass
    def make_var(self,var):
        self.init_varwithdims(var,[])
                                
# general class for reduced statistics based on 3d output
class statgroup_reduced(ncobject):    
    def __init__(self,outfile,description,logical):
        super(statgroup_reduced,self).__init__(outfile,description,logical)    
    def set_dims(self):
        self.init_vdims()
        self.init_hdims()
    def init_hdims(self):
        pass
    def make_var(self,var,mask=None):
        pass
        
class statgroupmean_xz(statgroup_reduced):    
    def __init__(self,outfile,description):
        super(statgroupmean_xz,self).__init__(outfile,description,ldiag_xz)    
    def init_hdims(self):
        self.init_xdims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,['z','x'],mask)
        
class statgroupmean_yz(statgroup_reduced):    
    def __init__(self,outfile,description):
        super(statgroupmean_yz,self).__init__(outfile,description,ldiag_yz)    
    def init_hdims(self):
        self.init_ydims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,['z','y'],mask)

class statgroupcross_xz(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_xz,self).__init__(outfile,description,lcross_xz)
        self.sel=sel    
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims(self.sel)
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])

class statgroupcross_yz(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_yz,self).__init__(outfile,description,lcross_yz)  
        self.sel=sel      
    def init_hdims(self):
        self.init_xdims(self.sel)
        self.init_ydims()
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])

class statgroupcross_xy(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_xy,self).__init__(outfile,description,lcross_xy)
        self.sel=sel      
    # overwrite the vertical to make a level selection 
    def set_dims(self):
        self.init_vdims(self.sel)
        self.init_hdims()
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims()
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])
      
class statgroupspectra(ncobject):        
    def __init__(self,outfile,description):
        super(statgroupspectra,self).__init__(outfile,description,lspec)
        self.initiated=False
    def put_make_var(self,var,field):
        if self.active:
            for specbounds in spectralevels:
                self.get_spacing()
                speclower=specbounds[0]
                specupper=specbounds[1]
                levelstring='_'+str(speclower)+'_'+str(specupper)
                # calculate distance between grid points
                if self.direction=='y':
                    nsegmentsmax=shape(field)[0]*(specupper-speclower)   
                    dataout=zeros((nsegmentsmax,shape(field)[1]))
                    prepare_spectra_y(dataout,field,shape(field)[0],shape(field)[1],speclower,specupper)
                elif self.direction=='x':
                    nsegmentsmax=shape(field)[1]*(specupper-speclower)   
                    dataout=zeros((nsegmentsmax,shape(field)[0]))
                    prepare_spectra_x(dataout,field,shape(field)[0],shape(field)[1],speclower,specupper)
                (p,wavenr)=areaspectra.spectrum_peri(dataout, Fs=1/self.dgrid, pad=False, smooth=False,rmzf=True,scale_by_freq=True)
                if(self.initiated==False):
                    self.init_dim('wavenr','wave number [m-1]',wavenr)
                    self.initiated=True
                if(self.tstep==0):
                    if 'ze' in allfields[var]:
                        self.ze=self.gdim('z')
                        self.init_var(var+levelstring,'power spectrum ('+self.direction+'-direction) of '+var+' at '+str(int(self.ze[speclower]))+'-'+str(int(self.ze[specupper]))+' meter (levels '+str(speclower)+'-'+str(specupper)+')','PSD',('time','wavenr'))
                    else:
                        zmin=self.gdim('z')[:-1]
                        zplus=self.gdim('z')[1:]
                        zhalf=0.5*(zmin+zplus)
                        bottom=-zhalf[0]
                        self.zc=hstack(([bottom],zhalf))
                        self.init_var(var+levelstring,'power spectrum ('+self.direction+'-direction) of '+var+' at '+str(int(self.zc[speclower]))+'-'+str(int(self.zc[specupper]))+' meter (levels '+str(speclower)+'-'+str(specupper)+')','PSD',('time','wavenr'))
                self.put_var(var+levelstring,p)
                
class statgroupspectra_y(statgroupspectra):        
    def get_spacing(self):
        y=self.gdim('y')
        self.dgrid=y[1]-y[0]
        self.direction='y'
        
class statgroupspectra_x(statgroupspectra):        
    def get_spacing(self):
        x=self.gdim('x')
        self.dgrid=x[1]-x[0]
        self.direction='x'

class statgroupclouds(statgroup_reduced):    
    def __init__(self,outfile,description):
        super(statgroupclouds,self).__init__(outfile,description,lclouds)
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims()
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])
    def init_var(self,var,longname,units,dims):
        if var in ['QV','QC','QR','QS','QG','QI']:
             so=self.outfile.createVariable(var, 'f4', dims,zlib=lzlib,least_significant_digit=6)
             so.missing_value = nan
             so.long_name=longname
             so.units=units
        else:
             so=self.outfile.createVariable(var, 'f4', dims,zlib=lzlib,least_significant_digit=3)
             so.missing_value = nan
             so.long_name=longname
             so.units=units
                             
# class for postprocessing data
class dataprocessor(get_variable_class):
    def __init__(self):
        self.data=[]
        self.helper=[]
        self.clouds=[]
        self.stat_1d=statgroup_1d('stat_1d.'+case+'.nc','MONC mean profile diagnostics')
        self.samp_1d=statgroup_1d('samp_1d.'+case+'.nc','MONC sampled profile diagnostics')
        self.cross_xz=statgroupcross_xz('cross_xz.'+case+'.nc','MONC cross-sections in the x,z plane',ysel)
        self.cross_yz=statgroupcross_yz('cross_yz.'+case+'.nc','MONC cross-sections in the y,z plane',xsel)
        self.cross_xy=statgroupcross_xy('cross_xy.'+case+'.nc','MONC cross-sections in the x,y plane',zsel)
        self.stat_xz=statgroupmean_xz('stat_xz.'+case+'.nc','MONC mean diagnostics in the x,z plane')
        self.stat_yz=statgroupmean_yz('stat_yz.'+case+'.nc','MONC mean diagnostics in the y,z plane')
        self.stat_int=statgroup_int('stat_int.'+case+'.nc','MONC column integrated diagnostics')
        self.stat_dom=statgroup_dom('stat_dom.'+case+'.nc','MONC domain integrated diagnostics')
        self.spec_x=statgroupspectra_x('spec_x.'+case+'.nc','MONC spectra along the x-direction')
        self.spec_y=statgroupspectra_y('spec_y.'+case+'.nc','MONC spectra along the y-direction')
        if lsamp:
            self.init_masks(['cld','cldupd','cldupdw1','upd'])
    def calc_masks(self):
        if lsamp:
            w=self.gv('w')
            cldze=self.helper.cld # should be updated with interpolated saturation deficit
            self.masks['cld'].setfieldze(cldze)
            self.masks['cldupd'].setfieldze(cldze*(w>0.0))
            self.masks['cldupdw1'].setfieldze(cldze*(w>1.0))
            self.masks['upd'].setfieldze((w>0.0))
            del cldze
            self.masks['cld'].setfield(self.helper.cld)
            self.masks['cldupd'].setfield(self.helper.cld*(self.helper.wzc>0.0))
            self.masks['cldupdw1'].setfield(self.helper.cld*(self.helper.wzc>1.0))
            self.masks['upd'].setfield((self.helper.wzc>0.0))
            cldxe=self.helper.cld # should be updated with interpolated saturation deficit
            wxe=0.5*(self.helper.wzc[:,:,:]+vstack((self.helper.wzc[-1,:,:][None,:,:],self.helper.wzc[:-1,:,:])))
            self.masks['cld'].setfieldxe(cldxe)
            self.masks['cldupd'].setfieldxe(cldxe*(wxe>0.0))
            self.masks['cldupdw1'].setfieldxe(cldxe*(wxe>1.0))
            self.masks['upd'].setfieldxe((wxe>0.0))
            del cldxe,wxe
            cldye=self.helper.cld # should be updated with interpolated saturation deficit
            wye=0.5*(self.helper.wzc[:,:,:]+hstack((self.helper.wzc[:,-1,:][:,None,:],self.helper.wzc[:,:-1,:])))
            self.masks['cld'].setfieldye(cldye)
            self.masks['cldupdw1'].setfieldye(cldye*(wye>1.0))
            self.masks['cldupd'].setfieldye(cldye*(wye>0.0))
            self.masks['upd'].setfieldye((wye>0.0))
    def init_masks(self,masks):
        self.masks={}
        for maskname in masks:
            self.masks[maskname]=mask()
    def app_tstep(self,data):
        self.data=data
        self.varkeys=self.data.variables.keys()
    def app_tstep(self,data,helper):
        self.data=data
        self.helper=helper
        self.varkeys=self.data.variables.keys()
        timestep=self.data.variables['timestep'][0]
        self.calc_masks()
        self.clouds=statgroupclouds("clouds/clouds."+case+".%05d.nc" %timestep,'3d in-cloud variable fields at %05d seconds' %timestep)
        # opening and closing enables us to check the data already while it is being processed
        self.stat_1d.opener(data)
        self.samp_1d.opener(data)
        self.cross_xz.opener(data)
        self.cross_yz.opener(data)
        self.cross_xy.opener(data)
        self.stat_xz.opener(data)
        self.stat_yz.opener(data)
        self.stat_int.opener(data)
        self.stat_dom.opener(data)
        self.clouds.opener(data)
        self.spec_x.opener(data)
        self.spec_y.opener(data)
        self.processor()
        self.stat_1d.closer()
        self.samp_1d.closer()
        self.cross_xz.closer()
        self.cross_yz.closer()
        self.cross_xy.closer()
        self.stat_xz.closer()
        self.stat_yz.closer()
        self.stat_int.closer()
        self.stat_dom.closer()
        self.clouds.closer()
        self.spec_x.closer()
        self.spec_y.closer()
    def stat_var(self,var,field):
        self.stat_1d.put_make_var(var,mean_2d(field))      
        self.stat_xz.put_make_var(var,mean_xz(field))
        self.stat_yz.put_make_var(var,mean_yz(field))
        self.masked_var(var,field)
    def process_var(self,var,field):
        self.stat_1d.put_make_var(var,mean_2d(field))
        self.cross_xz.put_make_var(var,field[:,ysel,:])
        self.cross_yz.put_make_var(var,field[xsel,:,:]) 
        self.cross_xy.put_make_var(var,field[:,:,zsel])           
        self.stat_xz.put_make_var(var,mean_xz(field))
        self.stat_yz.put_make_var(var,mean_yz(field))
        if 'makespectra' in allfields[var]:
            self.spec_x.put_make_var(var,field)
            self.spec_y.put_make_var(var,field)
        if var in progfields:
            temp=1.0*field #force copy
            whereinf=(self.helper.cld==0);
            temp[whereinf] = nan
            self.clouds.put_make_var(var,temp)
        self.masked_var(var,field)
    def ref_var(self,var,field):
        self.stat_1d.put_make_var(var,field)
    def masked_var(self,var,field):
        for mask in self.masks.keys():         
            if 'xe' in allfields[var]:
                self.samp_1d.put_make_sampvar(var,mean_2d(field*self.masks[mask].fieldxe)/mean_2d(self.masks[mask].fieldxe),mask)           
            elif 'ye' in allfields[var]:
                self.samp_1d.put_make_sampvar(var,mean_2d(field*self.masks[mask].fieldye)/mean_2d(self.masks[mask].fieldye),mask)  
            elif 'ze' in allfields[var]:
                self.samp_1d.put_make_sampvar(var,mean_2d(field*self.masks[mask].fieldze)/mean_2d(self.masks[mask].fieldze),mask)                      
            else:
                self.samp_1d.put_make_sampvar(var,mean_2d(field*self.masks[mask].field)/mean_2d(self.masks[mask].field),mask)
    def int_var(self,var,field):
        self.stat_int.put_make_var(var,field)
        if 'max' in allfields[var]:
           self.stat_dom.put_make_var(var,nanmax(field))
        elif 'min' in allfields[var]:
           self.stat_dom.put_make_var(var,nanmin(field))
        else:
           self.stat_dom.put_make_var(var,mean_2d(field))
    def dom_var(self,var,value):
        self.stat_dom.put_make_var(var,value)
    # ACTUAL CALCULATIONS HAPPEN HERE
    def processor(self):
        u=self.gv('u')
        v=self.gv('v')
        w=self.gv('w')
        deltheta=self.gv('th')
        thetaref=self.gref('thref')
        delp=self.gv('p')
        qv=self.gq('q',0)
        qc=self.gq('q',1)
        pref=self.gref('prefn')
        p=delp+pref[None,None,:]
        exn=(pref/psfr)**(rd/cp)
        theta=thetaref[None,None,:]+deltheta
        t=theta*exn[None,None,:]
        self.process_var('U',u)      
        self.process_var('V',v)      
        self.process_var('W',w)
        self.process_var('QC',qc)
        self.process_var('QV',qv)
        self.process_var('QT',qc+qv)
        self.process_var('THETA',theta)  
        self.process_var('P',p)
        self.process_var('T',t)
        self.process_var('TMSE',t+(grav/cp)*self.helper.zc[None,None,:]+(rlvap/cp)*qv)
        self.process_var('TLISE',t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*qc)
        rhon=self.gref('rhon')       
        rho=self.gref('rho')
        self.ref_var('EXNREF',exn)
        self.ref_var('RHOREF',rhon)
        self.ref_var('THETAREF',thetaref)
        self.ref_var('RHOREFH',rho)
        thetarhox=theta*(1+(rvord-1)*qv-qc)
        self.process_var('THETARHOX',thetarhox)
        buoyx=grav*deviation_2d(thetarhox)/thetaref[None,None,:]
        self.process_var('BUOYX',buoyx)  
        dbuoyx=deviation_2d(buoyx)
        self.stat_var('BUOYXVAR',dbuoyx*dbuoyx)      
        del dbuoyx
        del thetarhox,buoyx
        thetal=theta-(rlvap/(cp*exn))*qc
        self.process_var('THETAL',thetal) 
        del thetal
        qsat=qsa1/(p*exp(qsa2*(t-tk0c)/(t-qsa3))-qsa4) 
        qsati=qis1/(p*exp(qis2*(t-tk0c)/(t-qis3))-qis4) 
        self.process_var('QSAT',qsat) 
        self.process_var('QSATI',qsati) 
        self.process_var('RH',(qc+qv)/qsat)
        self.process_var('RHI',(qc+qv)/qsati) 
        del qsat,qsati
        # variances only produce statistics, not cross-sections     
        du=deviation_2d(u)
        self.stat_var('UVAR',du*du)      
        dv=deviation_2d(v)
        self.stat_var('VVAR',dv*dv)      
        dw=deviation_2d(w)
        self.stat_var('WVAR',dw*dw)      
        self.stat_var('TKERES',du*du+dv*dv+dw*dw)
        del du,dv,dw
        dth=deviation_2d(theta)
        self.stat_var('THETAVAR',dth*dth)    
        del dth
        dqv=deviation_2d(qv)
        self.stat_var('QVVAR',dqv*dqv)    
        del dqv
        dqt=deviation_2d(qv+qc)
        self.stat_var('QTVAR',dqt*dqt)
        del dqt    
        dp=deviation_2d(p)
        self.process_var('DP',dp)
        del dp       
        # height integrated variables 
        self.int_var('WMIN',self.helper.wmin)      
        self.int_var('WMAX',self.helper.wmax)
        self.int_var('CLDTOP',nanmax(self.helper.cld*self.clouds.zc[None,None,:],axis=2))   
        self.int_var('CLDW1TOP',nanmax(self.helper.cld*(self.helper.wzc>1.0)*self.clouds.zc[None,None,:],axis=2))
        # domain integrated variables 
        self.dom_var('CC',mean_2d(nanmax(self.helper.cld,axis=2)))      
        for mask in self.masks.keys():
            self.samp_1d.put_make_sampvar('frac',mean_2d(self.masks[mask].field),mask)                 
            self.samp_1d.put_make_sampvar('fracxe',mean_2d(self.masks[mask].fieldxe),mask)           
            self.samp_1d.put_make_sampvar('fracye',mean_2d(self.masks[mask].fieldye),mask)  
            self.samp_1d.put_make_sampvar('fracze',mean_2d(self.masks[mask].fieldze),mask)  
                                
# process 3d output fields
def process_3doutput():
    global heighthelper
    helper=nchelper()
    processor=dataprocessor()
    for file_to_process in filelist['3ddump']:
        moncdata=Dataset(file_to_process,'r',format='NETCDF4')
        helper.update(moncdata)
        processor.app_tstep(moncdata,helper)
        print 'time is '+str(time.clock()-start)
     
# replace missing values for reading into ncview
# and copy to project (storage) directory
def copy_files_to_project():
    scratchfiles=glob.glob(scratchdir+'*'+case+'*.nc')
    scratchfilescloud=glob.glob(scratchdir+'clouds/clouds.'+case+'*.nc')
    for scratchfile in scratchfiles+scratchfilescloud:
        print scratchfile
        scratchfiledata=Dataset(scratchfile,'r+',format='NETCDF4')
        for var in scratchfiledata.variables:
            if(type(scratchfiledata.variables[(var)])==numpy.ma.masked_array):
                vardata=scratchfiledata.variables[(var)][:]
                vardata=vardata.filled(nan)
                whereinf=isinf(vardata);
                vardata[whereinf]=nan
        scratchfiledata.close()
        print 'copying, time is '+str(time.clock()-start)
    for scratchfile in scratchfiles:        
        shutil.copy(scratchfile,projectdir)
    make_tarfile(scratchdir+'clouds.'+case+'.tar',glob.glob(scratchdir+'clouds/clouds.'+case+'*.nc'))
    shutil.copy(scratchdir+'clouds.'+case+'.tar',projectdir)

# update the variables to post-process by level type
def update_variables():
    global progfields,prevfields,derivedfields,intfields,vertfields,allfields   
    allfields={}
    for i in [progfields,prevfields,derivedfields,intfields,vertfields,domfields,fracfields]:
        allfields.update(i)
    alls=allfields.keys()

########### MAIN PROGRAM ########### 
def runme():
    update_variables()
    mkdir_p(scratchdir)
    mkdir_p(scratchdir+'/clouds')
    mkdir_p(projectdir)
    make_filelist()
    process_3doutput()
    copy_files_to_project()

# ACTUALLY CALLS THE SCRIPT FROM THE COMMAND LINE
# Using the construction with
# if __name__ == "__main__"
# makes sure we can import the separate routines
 
if __name__ == "__main__":
    global case,exper
    case=sys.argv[1]
    exper=sys.argv[2]
    fulldir=fullbase+case+'/output/'
    scratchdir=scratchbase+case+'/'
    projectdir=projectbase+case+'/'
    runme()
