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
                   
## IMPORTS

from numpy import *
from netCDF4 import Dataset
from scipy import weave # to embed code in C
from scipy.weave import converters
from fieldslist import * # a separate file contains the actual variable list
from contextlib import closing
from pp_MONC_constants import *
from pp_MONC_settings import *
from scipy.interpolate import splrep,spalde

import sys # system functions
import glob # a libary for regular expressions in file names
import shutil # for copying the files to destination
import areaspectra # a separate library for computing spectra
import tarfile # for compressing the cloud field data
import errno
import time
import numpy
   
start=time.clock() 
numpy.seterr(invalid='ignore') # don't wine about nans

case=''
exper=''
fulldir=''
scratchdir=''
projectdir=''

##  EXTRA DECLARATIONS FOR FAST SATURATION PHYSICS

numpsatarr=zeros(350,double) # saturation pressure values numerator array
denpsatarr=zeros(350,double) # saturation pressure values denominator array
psatarr=zeros(350,double) # saturation pressure values denominator array
psatarrabc=zeros((3,350),double) # saturation pressure spline array

for i in range(350):
    numpsatarr[i]=qsa1/(0.01*exp(qsa2*(i-tk0c)/(i-qsa3)))  
    denpsatarr[i]=qsa4/(0.01*exp(qsa2*(i-tk0c)/(i-qsa3)))
    psatarr[i]=0.01*exp(qsa2*(i-tk0c)/(i-qsa3))

psatspline = splrep(range(150,350),psatarr[150:350], k=2)
splinecoeffs=spalde(range(150,350),psatspline)
for i in range(150,350):
   psatarrabc[0,i]=splinecoeffs[i-150][0]
   psatarrabc[1,i]=splinecoeffs[i-150][1]
   psatarrabc[2,i]=splinecoeffs[i-150][2]

## FUNCTIONS
        
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
        filelist[i].sort(key=len)       
 
# create a tarfile (for cloud fields)
def make_tarfile(output_filename, infiles):
    with closing(tarfile.open(output_filename,'w')) as tar:
        for infile in infiles:
            tar.add(infile)
                          
                                            
# prepare spectral data into 2d arrays
# containing stretches along all points 
# that match sampling criteria
# this is for y-direction geometry
# uses C-code (weave) for speed
def prepare_spectra_y(dataout,data,imax,jmax,kmin,kmax,force=0):
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
    weave.inline(code,['dataout','data','imax','jmax','kmin','kmax'],type_converters=converters.blitz,compiler='gcc',force=force)

# The spectra with along x-direction
def prepare_spectra_x(dataout,data,imax,jmax,kmin,kmax,force=0):
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
    weave.inline(code,['dataout','data','imax','jmax','kmin','kmax'],type_converters=converters.blitz,compiler='gcc',force=force)
         
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

def interpolate_xe(field):
    fieldxe=0.5*(field[:,:,:]+vstack((field[-1,:,:][None,:,:],field[:-1,:,:])))
    return fieldxe

def interpolate_ye(field):
    fieldye=0.5*(field[:,:,:]+hstack((field[:,-1,:][:,None,:],field[:,:-1,:])))
    return fieldye

def interpolate_ze(field):
    fieldze=dstack((0.5*(field[:,:,:-1]+field[:,:,1:]),field[:,:,-1]+(field[:,:,-1]-field[:,:,-2])))
    return fieldze

def interpolate_ze_1d(array):
    arrayze=hstack((0.5*(array[:-1]+array[1:]),[array[-1]+(array[-1]-array[-2])]))
    return arrayze
                
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

# process 3d output fields
def process_3doutput(processor):
    global heighthelper
    helper=nchelper()
    print 'performing data check'
    for file_to_process in filelist['3ddump']:
        print file_to_process
        moncdata=Dataset(file_to_process,'r',format='NETCDF4')
        moncdata.close()
    for file_to_process in filelist['3ddump']:
        print 'processing '+file_to_process
        moncdata=Dataset(file_to_process,'r',format='NETCDF4')
        helper.update(moncdata)
        processor.app_tstep(moncdata,helper)
        print 'cpu time is '+str(time.clock()-start)
        moncdata.close()

# replace missing values for reading into ncview
# and copy to project (storage) directory
def copy_files_to_project():
    scratchfiles=glob.glob(scratchdir+'*'+case+'*.nc')
    scratchfilescloud=glob.glob(scratchdir+'clouds/clouds.'+case+'*.nc')
    for scratchfile in scratchfiles+scratchfilescloud:
        print 'revising '+scratchfile
        print 'cpu time is '+str(time.clock()-start)
        scratchfiledata=Dataset(scratchfile,'r+',format='NETCDF4')
        for var in scratchfiledata.variables:
            if(type(scratchfiledata.variables[(var)])==numpy.ma.masked_array):
                vardata=scratchfiledata.variables[(var)][:]
                vardata=vardata.filled(nan)
                whereinf=isinf(vardata);
                vardata[whereinf]=nan
        scratchfiledata.close()
    for scratchfile in scratchfiles:
        print 'copying scratch file '+scratchfile+' to '+projectdir
        print 'cpu time is '+str(time.clock()-start)
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

def fastmoistphysics(tlisein,qtin,z,p,force=0):
    support = """
    #include <math.h>
    #include <stdio.h>
    """
    imax=shape(tlisein)[0]
    jmax=shape(tlisein)[1]
    kmax=shape(tlisein)[2]
    qc=zeros((imax,jmax,kmax),float32)
    tv=zeros((imax,jmax,kmax),float32)
    code = """
    int const imax="""+str(imax)+""";
    int const jmax="""+str(jmax)+""";
    int const kmax="""+str(kmax)+""";
    double const cp="""+str(cp)+""";
    double const rlvap="""+str(rlvap)+""";
    double const grav="""+str(grav)+""";
    double const rvord="""+str(rvord)+""";
    double const invrlvap=1.0/rlvap;
    double const invcp=1.0/cp;
    double const qsa1="""+str(qsa1)+""";
    double const qsa4="""+str(qsa4)+""";
    
    int niter;
    bool lunsat;
    double zt, ztold, zttry, zqt, ztlise, zz, zp;
    double znumpint, zdenpint, zpsat, zqsat;
    double zqc, ztv;
    double ztliseguess,ztliseguesstry;
    int ztint;
    
    for (int i=0; i<imax; ++i) {
    for (int j=0; j<jmax; ++j) {
    for (int k=1; k<kmax; ++k) {
      ztlise=tlisein(i,j,k);
      zqt=qtin(i,j,k);
      zz=z(k);
      zp=p(k);
      
      // start loop
      lunsat=true;
            
      //first temperature guess: dry parcel
      zt=ztlise-(grav*invcp)*zz;
      ztint=std::min(int(zt),349);
      znumpint=numpsatarr(ztint);
      zdenpint=denpsatarr(ztint);
      
      if(zqt*(zp-zdenpint)>znumpint) {
        lunsat=false; // definitely unsaturated
      }
      if(lunsat==false) {
        zpsat=psatarrabc(0,ztint)+(zt-ztint)*psatarrabc(1,ztint)+0.5*(zt-ztint)*(zt-ztint)*psatarrabc(2,ztint);
        if(not(zqt*(zp*zpsat-qsa4)>qsa1)) {
          lunsat=true; // definitely unsaturated
        }        
      }
     
      // the long saturation loop is needed: jaiks
      if(lunsat==false) {
       
         ztold=zt;
         niter=0;
         zqsat=qsa1/(zp*zpsat-qsa4);
         ztliseguess=zt+(grav*invcp)*zz-(rlvap*invcp)*std::max(zqt-zqsat,0.0);
                  
         // Calculate Newton-Raphson iteration at perturbed temperature
         // Follow the saturation curves at the previous temperature
         zttry=zt-0.002;

         zpsat=zpsat-0.002*psatarrabc(1,ztint)+0.5*4.0e-6*psatarrabc(2,ztint);
         zqsat=qsa1/(zp*zpsat-qsa4);
         ztliseguesstry =zttry+(grav*invcp)*zz-(rlvap*invcp)*std::max(zqt-zqsat,0.0);
         
         zt = std::min(std::max(ztlise-(grav*invcp)*zz,zt-(ztliseguess-ztlise)/((ztliseguess-ztliseguesstry)*500.0)),ztlise-(grav*invcp)*zz+(rlvap*invcp)*zqt);
                  
         while((std::abs(zt-ztold) > 0.002) and (niter<40)) {
           niter = niter+1;
           ztold=zt;
           ztint=std::min(int(zt),349);
           zpsat=psatarrabc(0,ztint)+(zt-ztint)*psatarrabc(1,ztint)+0.5*(zt-ztint)*(zt-ztint)*psatarrabc(2,ztint);
           zqsat=qsa1/(zp*zpsat-qsa4);
           ztliseguess=zt+(grav*invcp)*zz-(rlvap*invcp)*std::max(zqt-zqsat,0.0);
                    
           // Calculate Newton-Raphson iteration at perturbed temperature
           // Follow the saturation curves at the previous temperature
           zttry=zt-0.002;
           zpsat=zpsat-0.002*psatarrabc(1,ztint)+0.5*4.0e-6*psatarrabc(2,ztint);
           zqsat=qsa1/(zp*zpsat-qsa4);
           ztliseguesstry =zttry+(grav*invcp)*zz-(rlvap*invcp)*std::max(zqt-zqsat,0.0);                

           zt = std::min(std::max(ztlise-(grav*invcp)*zz,zt-(ztliseguess-ztlise)/((ztliseguess-ztliseguesstry)*500.0)),ztlise-(grav*invcp)*zz+(rlvap*invcp)*zqt);
        }
                 
        zqc=cp*invrlvap*(zt-ztlise+(grav*invcp)*zz);
        
        if(zqc>0.0) {
          // moist parcel, adjust qc
          // conservation of lise, qt
          ztv=zt*(1.0+(rvord-1.0)*(zqt-zqc)-zqc);
        }        
        else {
          // parcel is dry after all, see below
          lunsat=true;
        }       
      }
      
      if(lunsat==true) {
        ztv=zt*(1.0+(rvord-1.0)*zqt);
        zqc=0.0;
      }
      
      tv(i,j,k)=ztv;
      qc(i,j,k)=zqc;
    }
    }
    }

    for (int i=0; i<imax; ++i) {
    for (int j=0; j<jmax; ++j) {
      tv(i,j,0)=tv(i,j,1);
      qc(i,j,0)=0.0;
    }
    }
    """
    weave.inline(code, ['tlisein','qtin','p','z','tv','qc','numpsatarr','denpsatarr','psatarrabc'],type_converters = converters.blitz,support_code=support,compiler='gcc',force=force)
    return tv,qc
    
## CLASSES

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
            return(self.data.variables[('pref')][:]*nan) 
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
                if index<len(self.data.variables[(key)]):
                    return self.data.variables[(key)][index,:,:,:]
                else:
                    return zeros(shape(self.data.variables[(key)][0]))
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

# a class for netcdf output                              
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
        self.force=1
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
                    prepare_spectra_y(dataout,field,shape(field)[0],shape(field)[1],speclower,specupper,self.force)
                elif self.direction=='x':
                    nsegmentsmax=shape(field)[1]*(specupper-speclower)   
                    dataout=zeros((nsegmentsmax,shape(field)[0]))
                    prepare_spectra_x(dataout,field,shape(field)[0],shape(field)[1],speclower,specupper,self.force)
                (p,wavenr)=areaspectra.spectrum_peri(dataout, Fs=1/self.dgrid, pad=False, smooth=False,rmzf=True,scale_by_freq=True)
                self.force=0
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
class dataorganizer(get_variable_class):
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
            self.init_masks(['cld','cldupd','cldupdw1','upd','buoyx','cldbuoyx'])
        self.force=1
    def calc_masks(self):
        if lsamp:
            deltheta=self.gv('th')
            thetaref=self.gref('thref')     
            delp=self.gv('p')
            pref=self.gref('prefn')
            exn=(pref/psfr)**(rd/cp)
            theta=thetaref[None,None,:]+deltheta
            t=theta*exn[None,None,:]
            del theta,exn,deltheta,thetaref	
	    qc=self.gq('q',nqc)
	    qi=self.gq('q',nqi)
	    qv=self.gq('q',nqv)
            tlise=t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*self.gq('q',nqc)
            qt=qc+qi+qv
	    tv=t*(1+(rvord-1)*qv-qc-qi)
            del qv,qc,qi,t
	    meantv=mean_2d(tv)
	    dtv=tv-meantv[None,None,:]
	    #tvcheck,qccheck=fastmoistphysics(tlise,qt,self.helper.zc,pref)
            #cldcheck=(qccheck>0.0)
            #dtvcheck=tvcheck-mean_2d(tv)
            #print mean_2d(dtvcheck)
            cldcheck=self.helper.cld
	    self.masks['cld'].setfield(cldcheck)
            self.masks['cldupd'].setfield(cldcheck*(self.helper.wzc>0.0))
            self.masks['cldupdw1'].setfield(cldcheck*(self.helper.wzc>1.0))
            self.masks['upd'].setfield((self.helper.wzc>0.0))
            self.masks['buoyx'].setfield((dtv>0.0))
            self.masks['cldbuoyx'].setfield((dtv>0.0)*cldcheck)   
	    del dtv,cldcheck
            # PHYSICS AT XE
            wxe=interpolate_xe(self.helper.wzc)
            tlisexe=interpolate_xe(tlise)
            qtxe=interpolate_xe(qt)
            tvxe,qcxe=fastmoistphysics(tlisexe,qtxe,self.helper.zc,pref,force=self.force)
	    self.force=0
	    del qtxe,tlisexe
	    dtvxe=tvxe-meantv
	    cldxe=(qcxe>1.0e-6)
            self.masks['cld'].setfieldxe(cldxe)
            self.masks['cldupd'].setfieldxe(cldxe*(wxe>0.0))
            self.masks['cldupdw1'].setfieldxe(cldxe*(wxe>1.0))
            self.masks['upd'].setfieldxe((wxe>0.0))
            self.masks['buoyx'].setfieldxe((dtvxe>0.0))
            self.masks['cldbuoyx'].setfieldxe((dtvxe>0.0)*cldxe)
            del cldxe,wxe,tvxe,dtvxe
	    # PHYSICS AT YE	    
            wye=interpolate_ye(self.helper.wzc)
            tliseye=interpolate_ye(tlise)
            qtye=interpolate_ye(qt)
            tvye,qcye=fastmoistphysics(tliseye,qtye,self.helper.zc,pref)
	    del qtye,tliseye
	    dtvye=tvye-meantv
	    cldye=(qcye>1.0e-6)
            self.masks['cld'].setfieldye(cldye)
            self.masks['cldupdw1'].setfieldye(cldye*(wye>1.0))
            self.masks['cldupd'].setfieldye(cldye*(wye>0.0))
            self.masks['upd'].setfieldye((wye>0.0))
            self.masks['buoyx'].setfieldye((dtvye>0.0))
            self.masks['cldbuoyx'].setfieldye((dtvye>0.0)*cldye)
            del cldye,wye,tvye,dtvye
            # PHYSICS AT ZE, TAKE INTO ACCOUNT SURFACE
            w=self.gv('w')
            qtze=interpolate_ze(qt)
	    tliseze=interpolate_ze(tlise)
	    prefze=interpolate_ze_1d(pref)
            tvze,qcze=fastmoistphysics(tliseze,qtze,self.gdim('z'),prefze)
	    del qtze,tliseze
	    dtvze=deviation_2d(tvze)
	    cldze=(qcze>1.0e-6)
            self.masks['cld'].setfieldze(cldze)
            self.masks['cldupd'].setfieldze(cldze*(w>0.0))
            self.masks['cldupdw1'].setfieldze(cldze*(w>1.0))
            self.masks['upd'].setfieldze((w>0.0))    
            self.masks['buoyx'].setfieldze((dtvze>0.0))
            self.masks['cldbuoyx'].setfieldze((dtvze>0.0)*cldze)
            del cldze,w,tvze,dtvze
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
        if not lsamp:
            return
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
                                
