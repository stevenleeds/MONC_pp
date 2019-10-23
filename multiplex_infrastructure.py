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
from multiplex_settings import *
from scipy.interpolate import splrep,spalde

import sys # system functions
import glob # a libary for regular expressions in file names
import shutil # for copying the files to destination
import areaspectra # a separate library for computing spectra
import tarfile # for compressing the cloud field data
import errno
import time
import numpy
import os
import datetime
import re

import numpy as np
   
start=time.clock() 
numpy.seterr(invalid='ignore') # don't wine about nans

outputconfig=opconfig()
sysconfig=sconfig()

##  EXTRA DECLARATIONS FOR FAST SATURATION PHYSICS

numpsatarr=zeros(350,double) # saturation pressure values numerator array
denpsatarr=zeros(350,double) # saturation pressure values denominator array
psatarr=zeros(350,double) # saturation pressure values denominator array
psatarrabc=zeros((3,350),double) # saturation pressure spline array

for i in range(350):
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

# write a log file
def writefinished(inputtext,outputfile):
    outfile = open(outputfile,'wb')
    outfile.write(inputtext)
    outfile.write("\n")
    outfile.write("This case was copied "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    outfile.write("\n")
    outfile.write("\n")
    outfile.close()
    
# make a filelist of the 3D output
def make_filelist():
    global filelist
    types = {'3ddump':'*.nc'} # file type list, currently contains only 3D output to be postprocessed
    filelist={}
    for i in types:
        # if overwrite option given, replace the processed files
        if sysconfig.overwrite:
            done_list=glob.glob(sysconfig.fulldir+types[i]+'_done')
            for file_to_process in done_list:
                os.rename(file_to_process,file_to_process.replace('_done','')) 
        filelist[i]=glob.glob(sysconfig.fulldir+types[i])
        filelist[i].sort() 
        filelist[i].sort(key=len)       

# make a filelist with one element. Forget about the done marking in this 
def make_onefile(filenumber):
    global filelist
    types = {'3ddump':'*.nc'} # file type list, currently contains only 3D output to be postprocessed
    filelist={}
    for i in types:
        filelist[i]=glob.glob(sysconfig.fulldir+types[i])
        filelist[i].sort() 
        filelist[i].sort(key=len)       
        filelist[i]=[filelist[i][filenumber]]
        print i,filelist[i]
 
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
    weave.inline(code,['dataout','data','imax','jmax','kmin','kmax'],type_convertfers=converters.blitz,compiler='gcc',force=force)

# mean value across 2d slab
def mean_xz(inputfield):
    meanfield=mean(inputfield,axis=1, dtype=numpy.float64)
    return meanfield
    
def mean_yz(inputfield):
    meanfield=mean(inputfield,axis=0, dtype=numpy.float64)
    return meanfield

def mean_2d(inputfield):
    meanfield1=mean(inputfield,axis=0,dtype=numpy.float64)
    meanfield=mean(meanfield1,axis=0,dtype=numpy.float64)
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

def xe_to_xc(field):
    fieldxc=0.5*(field[:,:,:]+vstack((field[:-1,:,:],field[-1,:,:][None,:,:])))    
    return fieldxc

def ye_to_yc(field):
    fieldyc=0.5*(field[:,:,:]+hstack((field[:,:-1,:],field[:,-1,:][:,None,:])))
    return fieldyc

def ze_to_zc(field):
    fieldzc=dstack((field[:,:,0]-0.5*(field[:,:,1]-field[:,:,0]),0.5*(field[:,:,:-1]+field[:,:,1:])))
    return fieldzc

def interpolate_xe(field):
    fieldxe=0.5*(field[:,:,:]+vstack((field[-1,:,:][None,:,:],field[:-1,:,:])))
    return fieldxe

def interpolate_ye(field):
    fieldye=0.5*(field[:,:,:]+hstack((field[:,-1,:][:,None,:],field[:,:-1,:])))
    return fieldye

def interpolate_ze(field):
    fieldze=dstack((0.5*(field[:,:,:-1]+field[:,:,1:]),field[:,:,-1]+0.5*(field[:,:,-1]-field[:,:,-2])))
    return fieldze

def interpolate_ze_1d(array):
    arrayze=hstack((0.5*(array[:-1]+array[1:]),[array[-1]+0.5*(array[-1]-array[-2])]))
    return arrayze

def cropped(array):
    if(outputconfig.nboundlines>0):
        return array[outputconfig.nboundlines:-outputconfig.nboundlines]
    else:
        return array
    
def dz(field,zin):  
    nanxy=float('nan')*field[:,:,0][:,:,None]
    return (concatenate((nanxy,(field[:,:,2:]-field[:,:,:-2])/(zin[None,None,2:]-zin[None,None,:-2]),nanxy),axis=2))
            
def dz_1d(field,zin):  
    nan_1d=float('nan')
    return (concatenate(([nan_1d],(field[2:]-field[:-2])/(zin[2:]-zin[:-2]),[nan_1d])))

# command to get a single variable from a file
def var_from_file(dataset,key):
    if(outputconfig.nboundlines>0):
        try:
            return dataset.variables[(key)][outputconfig.nboundlines:-outputconfig.nboundlines,outputconfig.nboundlines:-outputconfig.nboundlines,:]
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
        processor.app_tstep(moncdata,helper)
        print 'cpu time is '+str(time.clock()-start)
        moncdata.close()
        os.rename(file_to_process,file_to_process+'_done') 

# process 3d output fields
def process_onefile(processor):
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
        processor.app_tstep(moncdata,helper)
        print 'cpu time is '+str(time.clock()-start)
        moncdata.close()

# replace missing values for reading into ncview
# and copy to project (storage) directory
def copy_files_to_project():
    scratchfiles=glob.glob(sysconfig.scratchdir+'*'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber)
    scratchfilescloud=glob.glob(sysconfig.scratchdir+'pdfdata/pdfdata.*.%03d.nc'%sysconfig.filenumber)
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
    for scratchfile in scratchfiles+scratchfilescloud:
        print 'copying scratch file '+scratchfile+' to '+sysconfig.projectdir
        print 'cpu time is '+str(time.clock()-start)
        shutil.copy(scratchfile,sysconfig.projectdir)
    #make_tarfile(sysconfig.scratchdir+'clouds.'+sysconfig.exper+'.tar',glob.glob(sysconfig.scratchdir+'clouds/clouds.'+sysconfig.exper+'*.nc'))
    #shutil.copy(sysconfig.scratchdir+'clouds.'+sysconfig.exper+'.tar',sysconfig.projectdir)

# update the variables to post-process by level type
def update_variables():
    global progfields,prevfields,derivedfields,intfields,vertfields,allfields   
    allfields={}
    for i in [progfields,prevfields,derivedfields,intfields,vertfields,domfields,fracfields]:
        allfields.update(i)

def fastmoistphysics(tlisein,qtin,z,p,force=0):
    iter_count=0
    solved=False
    while((iter_count<20) and (solved==False)):
        try: 
            tv,qc=fastmoistphysics_inner(tlisein,qtin,z,p,force=0)
            solved=True
            print 'Succeeded at iteration '+str(iter_count)
        except:
            time.sleep(5)
            print 'Fastmoistphysics failed at iteration '+str(iter_count)
            if(iter_count==19):
                print 'Failed to do fastmoistphysics after 19 attempts' 
            iter_count=iter_count+1
            solved=False
    return tv,qc

def fastmoistphysics_inner(tlisein,qtin,z,p,force=0):
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
    double zpsat, zqsat;
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
      lunsat=false;
            
      //first temperature guess: dry parcel
      zt=ztlise-(grav*invcp)*zz;
      ztint=std::min(int(zt),349);

      zpsat=psatarrabc(0,ztint)+(zt-ztint)*psatarrabc(1,ztint)+0.5*(zt-ztint)*(zt-ztint)*psatarrabc(2,ztint);

      if(not(zqt*(zp*zpsat-qsa4)>qsa1)) {
        lunsat=true; // definitely unsaturated
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
    weave.inline(code, ['tlisein','qtin','p','z','tv','qc','psatarrabc'],type_converters = converters.blitz,support_code=support,compiler='gcc',force=force,verbose=2)
    return tv,qc
    
## CLASSES

# class for masks/conditional sampling
class mask:
    def __init__(self):
        pass
    def setfield(self,field):
        self.field=field
        self.ifrac=1.0/mean_2d(field)
    def setfieldxe(self,field):
        self.fieldxe=field
        self.ifracxe=1.0/mean_2d(field)
    def setfieldye(self,field):
        self.fieldye=field
        self.ifracye=1.0/mean_2d(field)
    def setfieldze(self,field):
        self.fieldze=field
        self.ifracze=1.0/mean_2d(field)
                
# a class to include some methods needed for both the help-variable class (which store temp variables) and the netcdf class
class get_variable_class():
    # gv being "get variable"
    def gv(self,key):
        if(outputconfig.nboundlines>0):
            if key in self.varkeys:
                return self.data.variables[(key)][self.step,outputconfig.nboundlines:-outputconfig.nboundlines,outputconfig.nboundlines:-outputconfig.nboundlines,:]
            else:
                return(self.data.variables[('p')][self.step,outputconfig.nboundlines:-outputconfig.nboundlines,outputconfig.nboundlines:-outputconfig.nboundlines,:]*nan)
        else:
            if key in self.varkeys:
                return self.data.variables[(key)][self.step,:,:,:]
            else:
                return(self.data.variables[('p')][self.step,:,:,:]*nan)          
    # gref being "get reference variable"
    def gref(self,key):
        if key in self.varkeys:
            return self.data.variables[(key)][self.step,:]
        else:
            return(self.data.variables[('pref')][:]*nan) 
    # gq "get moisture variable"
    def gq(self,key,index):
        if(outputconfig.nboundlines>0):
            if key in self.varkeys:
                if index<len(self.data.variables[(key)]) and index>-1:
                    return self.data.variables[(key)][self.step,index,outputconfig.nboundlines:-outputconfig.nboundlines,outputconfig.nboundlines:-outputconfig.nboundlines,:]
                else:
                    return zeros(shape(self.data.variables[(key)][self.step,0]))
            else:
                return(self.data.variables[('p')][self.step,outputconfig.nboundlines:-outputconfig.nboundlines,outputconfig.nboundlines:-outputconfig.nboundlines,:]*nan)
        else:
            if key in self.varkeys:
                if index<len(self.data.variables[(key)]) and index>-1:
                    return self.data.variables[(key)][self.step,index,:,:,:]
                else:
                    return zeros(shape(self.data.variables[(key)][self.step,0]))
            else:
                return(self.data.variables[('p')][self.step,:,:,:]*nan)    
    # gdim being "get dimension"
    def gdim(self,key):
        try:
            return self.data.variables[(key)][:]
        except:
            try:
                return array(range(len(self.data.dimensions[(key)])))
            except:
                return([nan]) 
    # gdim being "get dimension"
    def gdimt(self,key):
        try:
            return self.data.variables[(key)][0,:]
        except:
            try:
                return array(range(len(self.data.dimensions[(key)][0,:])))
            except:
                return([nan])
                
# a class to store derived variables from output which are needed relatively often       
class nchelper(object,get_variable_class):
    def __init__(self):
        self.data=[]
        self.varkeys=[]
        self.svlist=[]        
    def update(self,data,step):
        self.data=data
        self.step=step
        self.varkeys=self.data.variables.keys()
        w=self.gv('w')
        whalf=0.5*(w[:,:,1:]+w[:,:,:-1])
        bottom=-whalf[:,:,0]
        self.wzc=dstack((bottom[:,:,None],whalf))
        self.wmin=nanmin(w,axis=2)
        self.wmax=nanmax(w,axis=2)
        qci=self.gv('q_cloud_liquid_mass') # try to include ice
        self.cld=(qci>1.0e-6)
        self.cloudycolumn=1.0*(sum(self.cld,axis=2)>0)
        zmin=self.gdimt('z')[:-1]
        zplus=self.gdimt('z')[1:]
        zhalf=0.5*(zmin+zplus)
        bottom=-zhalf[0]
        self.zc=hstack(([bottom],zhalf))
        self.ze=self.gdimt('z')
        self.rhon=self.gref('rhon')
        self.xe=self.gdim('x')*outputconfig.dx
        self.ye=self.gdim('y')*outputconfig.dy
        if self.svlist==[]:
           self.initsvlist()
    def initsvlist(self):
        try:
            self.svlist=range(shape(self.data.variables[('q')])[0])
        except:
            self.svlist=[]
        for qnumber in [nqv,nqc,nqi,nqs,nqg,nqh]:
            try:
                self.svlist.remove(qnumber)
            except:
                pass
        for scalar_number in self.svlist:
            variable_to_add={
            'SCALAR%03d'%scalar_number:['SCALAR%03d'%scalar_number,u'-'],
            }
            allfields.update(variable_to_add)

# a class for netcdf output                              
class ncobject(object,get_variable_class):
    # class for writing to netcdf
    def __init__(self,outfile,description,logical=True):
        self.active=logical
        if self.active:
            self.data=[]
            self.outvars={}
            self.description=description
            self.ncoutname=sysconfig.scratchdir+outfile
            if sysconfig.overwrite:
                try:
                    os.remove(self.ncoutname)
                except:
                    pass
            if os.path.exists(self.ncoutname):
                self.outfile=Dataset(self.ncoutname,'a',format='NETCDF4',zlib=outputconfig.lzlib)
                self.outfile.history = self.outfile.history+'. Added ' + time.ctime(time.time()) +' '+myusername
            else:
                self.outfile=Dataset(self.ncoutname,'w',format='NETCDF4',zlib=outputconfig.lzlib)
                self.outfile.description=self.description+' for '+sysconfig.case+' '+sysconfig.exper
                self.outfile.history = 'Created ' + time.ctime(time.time())
                self.outfile.source = 'Created by user ' + myusername
                self.outfile.createDimension('time',0)
                timevar=self.outfile.createVariable('time', 'f8', ('time',),zlib=outputconfig.lzlib)
                setattr(timevar,'longname','time [s]')
            self.outfile.close()    
    def make_dims(self):
        pass
    def make_var(self):
        pass
    # LAYOUT OF POSTPROCESSING A TIME STEP
    def opener(self,data):
        if self.active:
            self.data=data      
            self.outfile=Dataset(self.ncoutname,'a',format='NETCDF4',zlib=outputconfig.lzlib)
    def opener2(self,data,step):
        self.step=step
        if self.active:
            self.t=0
            lenn=len(self.outfile.variables['time'])
            if(lenn>0):
                newtime=self.data.variables['time_series_600'][self.step]
                for tt in range(lenn,0,-1):
                    time_to_check=self.outfile.variables['time'][tt-1]
                    print 'time_to_check '+str(time_to_check)
                    print 'newtime '+str(newtime)
                    if(newtime>time_to_check):
                        self.t=tt
                        break
                if(self.t==0):
                    os.remove(self.ncoutname)
                    self.outfile=Dataset(self.ncoutname,'w',format='NETCDF4',zlib=outputconfig.lzlib)
                    self.outfile.description=self.description+' for '+sysconfig.case+' '+sysconfig.exper
                    self.outfile.history = 'Created ' + time.ctime(time.time())
                    self.outfile.source = 'Created by user ' + myusername
                    self.outfile.createDimension('time',0)
                    timevar=self.outfile.createVariable('time', 'f8', ('time',),zlib=outputconfig.lzlib)
                    setattr(timevar,'longname','time [s]')
            if(self.t==0):
                self.set_dims()
            self.outfile.variables['time'][self.t]=self.data.variables['time_series_600'][self.step]
    def closer(self):
        if self.active:
            self.outfile.close()
    # functions to initialize dimensions in output
    def set_dims(self):
        pass
    def init_dim(self,dimname,longdimname,dimvalues,sel=None):
        if sel==None:
            try:
                self.outfile.createDimension(dimname,len(dimvalues))
                var=self.outfile.createVariable(dimname, 'f8', (dimname,),zlib=outputconfig.lzlib)
                var[:]=dimvalues
            except:
                pass
        else:
            try:
                self.outfile.createDimension(dimname,len(sel))
                var=self.outfile.createVariable(dimname, 'f8', (dimname,),zlib=outputconfig.lzlib)
                var[:]=dimvalues[sel]
            except:
                pass
        try:
            setattr(var,'longname',longdimname)
        except:
            pass
    def init_dimxc(self,sel=None):
        xe=self.gdim('x')*outputconfig.dx
        xmin=hstack(([0],xe[:-1]))
        self.xc=cropped(0.5*(xmin+xe))
        self.init_dim('xc','x [m]',self.xc,sel)    
    def init_dimyc(self,sel=None):
        ye=self.gdim('y')*outputconfig.dy
        ymin=hstack(([0],ye[:-1]))
        self.yc=cropped(0.5*(ymin+ye))
        self.init_dim('yc','y [m]',self.yc,sel)
    def init_dimzc(self,sel=None):                      
        zmin=self.gdimt('z')[:-1]
        zplus=self.gdimt('z')[1:]
        zhalf=0.5*(zmin+zplus)
        bottom=-zhalf[0]
        self.zc=hstack(([bottom],zhalf))
        self.init_dim('zc','height [m]',self.zc,sel)    
    def init_dimxe(self,sel=None):
        self.xe=cropped(self.gdim('x')*outputconfig.dx)  
        self.init_dim('xe','x [m] (staggered)',self.xe,sel)
    def init_dimye(self,sel=None):
        self.ye=cropped(self.gdim('y')*outputconfig.dy)
        self.init_dim('ye','y [m] (staggered)',self.ye,sel)
    def init_dimze(self,sel=None):
        self.ze=self.gdimt('z')
        self.init_dim('ze','height [m] (bottom staggered)',self.ze,sel)
    def init_varwithdims(self,var,dimsin,mask=None):
        # initialise a variable with correct dimensions on a staggered grid
        dimsarr=['time']
        scalregex = re.compile('SCALAR...')
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
        so=self.outfile.createVariable(var, 'f4', dims,zlib=outputconfig.lzlib)
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
            if(self.t==0):
                self.make_var(var)
            self.put_var(var,field)
    def put_make_sampvar(self,var,field,mask):
        if self.active:
            if(self.t==0):
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
        super(statgroup_int,self).__init__(outfile,description,outputconfig.ldiag_int)    
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
        super(statgroup_dom,self).__init__(outfile,description,outputconfig.ldiag_int)    
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
        super(statgroupmean_xz,self).__init__(outfile,description,outputconfig.ldiag_xz)    
    def init_hdims(self):
        self.init_xdims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,['z','x'],mask)
        
class statgroupmean_yz(statgroup_reduced):    
    def __init__(self,outfile,description):
        super(statgroupmean_yz,self).__init__(outfile,description,outputconfig.ldiag_yz)    
    def init_hdims(self):
        self.init_ydims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,['z','y'],mask)

class statgroupcross_xz(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_xz,self).__init__(outfile,description,outputconfig.lcross_xz)
        self.sel=sel    
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims(self.sel)
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])

class statgroupcross_yz(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_yz,self).__init__(outfile,description,outputconfig.lcross_yz)  
        self.sel=sel      
    def init_hdims(self):
        self.init_xdims(self.sel)
        self.init_ydims()
    def make_var(self,var):
        self.init_varwithdims(var,['z','y','x'])

class statgroupcross_xy(statgroup_reduced):    
    def __init__(self,outfile,description,sel):
        super(statgroupcross_xy,self).__init__(outfile,description,outputconfig.lcross_xy)
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
        super(statgroupspectra,self).__init__(outfile,description,outputconfig.lspec)
        self.initiated=False
        self.force=1
    def put_make_var(self,var,field):
        if self.active:
            for speclevel in range(len(outputconfig.spectralevelsbot)):
                self.get_spacing()
                speclower=outputconfig.spectralevelsbot[speclevel]
                specupper=outputconfig.spectralevelstop[speclevel]
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
                if(self.t==0):
                    if 'ze' in allfields[var]:
                        self.ze=self.gdimt('z')
                        self.init_var(var+levelstring,'power spectrum ('+self.direction+'-direction) of '+var+' at '+str(int(self.ze[speclower]))+'-'+str(int(self.ze[specupper]))+' meter (levels '+str(speclower)+'-'+str(specupper)+')','PSD',('time','wavenr'))
                    else:
                        zmin=self.gdimt('z')[:-1]
                        zplus=self.gdimt('z')[1:]
                        zhalf=0.5*(zmin+zplus)
                        bottom=-zhalf[0]
                        self.zc=hstack(([bottom],zhalf))
                        self.init_var(var+levelstring,'power spectrum ('+self.direction+'-direction) of '+var+' at '+str(int(self.zc[speclower]))+'-'+str(int(self.zc[specupper]))+' meter (levels '+str(speclower)+'-'+str(specupper)+')','PSD',('time','wavenr'))
                self.put_var(var+levelstring,p)
                
class statgroupspectra_y(statgroupspectra):        
    def get_spacing(self):
        y=self.gdim('y')*outputconfig.dy
        self.dgrid=y[1]-y[0]
        self.direction='y'
        
class statgroupspectra_x(statgroupspectra):        
    def get_spacing(self):
        x=self.gdim('x')*outputconfig.dx
        self.dgrid=x[1]-x[0]
        self.direction='x'

class statgroupclouds(statgroup_reduced):    
    def __init__(self,outfile,description):
        super(statgroupclouds,self).__init__(outfile,description,outputconfig.lclouds)
    def opener(self,data):
        super(statgroupclouds,self).opener(data)
    def init_hdims(self):
        self.init_xdims()
        self.init_ydims()
    def make_var(self,var,mask=None):
        self.init_varwithdims(var,['z','y','x'],mask)
    def init_var(self,var,longname,units,dims):
        if any(text in var for text in ['QV','QC','QR','QS','QG','QI']):
            so=self.outfile.createVariable(var, 'f4', dims,zlib=outputconfig.lzlib,least_significant_digit=8)
            so.missing_value = nan
            so.long_name=longname
            so.units=units
        else:
            so=self.outfile.createVariable(var, 'f4', dims,zlib=outputconfig.lzlib,least_significant_digit=5)
            so.missing_value = nan
            so.long_name=longname
            so.units=units
                             
# class for postprocessing data
class dataorganizer(get_variable_class):
    def __init__(self):
        self.data=[]
        self.helper=[]
        self.clouds=[]
        self.stat_1d=statgroup_1d('stat_1d.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC mean profile diagnostics')
        self.samp_1d=statgroup_1d('samp_1d.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC sampled profile diagnostics')
        self.cross_xz=statgroupcross_xz('cross_xz.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC cross-sections in the x,z plane',outputconfig.ysel)
        self.cross_yz=statgroupcross_yz('cross_yz.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC cross-sections in the y,z plane',outputconfig.xsel)
        self.cross_xy=statgroupcross_xy('cross_xy.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC cross-sections in the x,y plane',outputconfig.zsel)
        self.stat_xz=statgroupmean_xz('stat_xz.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC mean diagnostics in the x,z plane')
        self.stat_yz=statgroupmean_yz('stat_yz.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC mean diagnostics in the y,z plane')
        self.stat_int=statgroup_int('stat_int.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC column integrated diagnostics')
        self.stat_dom=statgroup_dom('stat_dom.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC domain integrated diagnostics')
        self.spec_x=statgroupspectra_x('spec_x.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC spectra along the x-direction')
        self.spec_y=statgroupspectra_y('spec_y.'+sysconfig.exper+'.%03d.nc'%sysconfig.filenumber,'MONC spectra along the y-direction')
        if outputconfig.lsamp:
            self.init_masks(['tr_mean','tr_cc','tr_sh','tr_env','br_mean','br_cc','br_sh','br_env','bl_mean','bl_cc','bl_sh','bl_env','tl_mean','tl_cc','tl_sh','tl_env',])
        self.force=1
    def calc_masks(self):
        if outputconfig.lsamp:
            deltheta=self.gv('th')
            thetaref=self.gref('thref')   
            w=self.gv('w')  
            delp=self.gv('p')
            pref=self.gref('prefn')
            exn=(pref/psfr)**(rd/cp)
            theta=thetaref[None,None,:]+deltheta
            t=theta*exn[None,None,:]
            del theta,exn,thetaref
            qc=self.gv('q_cloud_liquid_mass')
            qi=0.0*qc
            qv=self.gv('q_vapour')
            tlise=t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*qc-(rlsub/cp)*qi
            qt=qc+qv+qi
            tv=t*(1+(rvord-1)*qv-qc-qi)
            meantv=mean_2d(tv)
            
            del qv,qi,t,deltheta
              
            quadrants = ['tr','tl','bl','br']
            segments = ['mean','cc','sh','env']

            for q in quadrants:
                # Define masks for each quadrant
                a = np.zeros((160, 160, 161))
                nx, ny, nz = np.shape(a)

                # Assign values to a
                if q == 'tr':
                    a[nx/2:nx, ny/2:ny, :] = 1
                elif q == 'tl':
                    a[0:nx/2, ny/2:ny, :] = 1
                elif q == 'bl':
                    a[0:nx/2, 0:ny/2, :] = 1
                elif q == 'br': 
                    a[nx/2:nx, 0:ny/2, :] = 1
                
                for s in segments:
                    dtv=tv-meantv[None,None,:]
                    cldcheck=self.helper.cld  
                    q_purity = self.gv('q_purity_tracer')

                    # Physics at xe
                    wxe=interpolate_xe(self.helper.wzc)
                    tlisexe=interpolate_xe(tlise)
                    qtxe=interpolate_xe(qt)
                    tvxe,qcxe=fastmoistphysics(tlisexe,qtxe,self.helper.zc,pref,force=self.force)
                    self.force=0
                    del qtxe,tlisexe
                    dtvxe=tvxe-meantv
                    cldxe=(qcxe>1.0e-6)
                    q_purityxe = interpolate_xe(q_purity)

                    # Physics at ye
                    wye=interpolate_ye(self.helper.wzc)
                    tliseye=interpolate_ye(tlise)
                    qtye=interpolate_ye(qt)
                    tvye,qcye=fastmoistphysics(tliseye,qtye,self.helper.zc,pref)
                    del qtye,tliseye
                    dtvye=tvye-meantv
                    cldye=(qcye>1.0e-6)
                    q_purityye = interpolate_ye(q_purity)

                    # Physics at ze
                    w=self.gv('w')
                    qtze=interpolate_ze(qt)
                    tliseze=interpolate_ze(tlise)
                    prefze=interpolate_ze_1d(pref)
                    tvze,qcze=fastmoistphysics(tliseze,qtze,self.gdimt('z'),prefze)
                    del qtze,tliseze
                    dtvze=deviation_2d(tvze)
                    cldze=(qcze>1.0e-6)
                    q_purityze = interpolate_ze(q_purity)

                    # Mean of quadrant
                    if s == 'mean':
                        self.masks[q+'_'+s].setfield(a > 0.5)
                        # Physics at xe
                        self.masks[q+'_'+s].setfieldxe(a > 0.5)
                        # Physics at ye
                        self.masks[q+'_'+s].setfieldye(a > 0.5)
                        # Physics at ze
                        self.masks[q+'_'+s].setfieldze(a > 0.5)

                    # Environment
                    if s == 'env':
                        self.masks[q+'_'+s].setfield(np.logical_and(a > 0.5, q_purity <= 1.0E-3))
                        del q_purity
                        # Physics at xe
                        self.masks[q+'_'+s].setfieldxe(np.logical_and(a > 0.5, q_purityxe <= 1.0E-3))
                        del cldxe,wxe,tvxe,dtvxe,q_purityxe
                        # Physics at ye
                        self.masks[q+'_'+s].setfieldye(np.logical_and(a > 0.5, q_purityye <= 1.0E-3))
                        del cldye,wye,tvye,dtvye,q_purityye
                        # Physics at ze
                        self.masks[q+'_'+s].setfieldze(np.logical_and(a > 0.5, q_purityze <= 1.0E-3))
                        del cldze,tvze,dtvze,q_purityze

                    # Cloud Core
                    elif s == 'cc':
                        attempt_1 = np.logical_and(a > 0.5, q_purity > 1.0E-3)
                        attempt_2 = np.logical_and(qc > 1.0E-5, w > 0.5)
                        self.masks[q+'_'+s].setfield(np.logical_and(attempt_1, attempt_2)) 
                        del attempt_1,attempt_2
                        # Physics at xe
                        attempt_1 = np.logical_and(a > 0.5, q_purityxe > 1.0E-3)
                        attempt_2 = np.logical_and(qcxe > 1.0E-5, wxe > 0.5)
                        self.masks[q+'_'+s].setfieldxe(np.logical_and(attempt_1, attempt_2))
                        del cldxe,wxe,tvxe,dtvxe,q_purityxe,qcxe,attempt_1,attempt_2
                        # Physics at ye
                        attempt_1 = np.logical_and(a > 0.5, q_purityye > 1.0E-3)
                        attempt_2 = np.logical_and(qcye > 1.0E-5, wye > 0.5)
                        self.masks[q+'_'+s].setfieldye(np.logical_and(attempt_1, attempt_2))
                        del cldye,wye,tvye,dtvye,q_purityye,qcye,attempt_1,attempt_2
                        # Physics at ze
                        attempt_1 = np.logical_and(a > 0.5, q_purityze > 1.0E-3)
                        attempt_2 = np.logical_and(qcze > 1.0E-5, w > 0.5)
                        self.masks[q+'_'+s].setfieldze(np.logical_and(attempt_1, attempt_2))
                        del cldze,tvze,dtvze,q_purityze,qcze,attempt_1,attempt_2,w
                    
                    # Shell
                    elif s == 'sh':
                        attempt_1 = np.logical_and(a > 0.5, q_purity*(self.helper.zc[None,None,:]>749.0) > 1.0E-3)
                        attempt_2 = np.logical_or(w <= 0.5, qc <=1.0E-5)
                        self.masks[q+'_'+s].setfield(np.logical_and(attempt_1, attempt_2))
                        del attempt_1,attempt_2
                        # Physics at xe
                        attempt_1 = np.logical_and(a > 0.5, q_purityxe*(self.helper.zc[None,None,:]>749.0) > 1.0E-3)
                        attempt_2 = np.logical_or(w <= 0.5, qcxe <=1.0E-5)
                        self.masks[q+'_'+s].setfieldxe(np.logical_and(attempt_1, attempt_2))
                        del cldxe,wxe,tvxe,dtvxe,q_purityxe,qcxe,attempt_1,attempt_2
                        # Physics at ye
                        attempt_1 = np.logical_and(a > 0.5, q_purityye*(self.helper.zc[None,None,:]>749.0) > 1.0E-3)
                        attempt_2 = np.logical_or(w <= 0.5, qcye <=1.0E-5)
                        self.masks[q+'_'+s].setfieldye(np.logical_and(attempt_1, attempt_2))
                        del cldye,wye,tvye,dtvye,q_purityye,qcye,attempt_1,attempt_2
                        # Physics at ze
                        attempt_1 = np.logical_and(a > 0.5, q_purityze*(self.helper.ze[None,None,:]>749.0) > 1.0E-3)
                        attempt_2 = np.logical_or(w <= 0.5, qcze <=1.0E-5)
                        self.masks[q+'_'+s].setfieldze(np.logical_and(attempt_1, attempt_2))
                        del cldze,tvze,dtvze,q_purityze,qcze,attempt_1,attempt_2,w

            del tv, dtv, qc

    def init_masks(self,masks):
        self.masks={}
        for maskname in masks:
            self.masks[maskname]=mask()
    def app_tstep(self,data,helper):
        self.data=data
        self.helper=helper
        self.varkeys=self.data.variables.keys()
        nsteps=len(self.data.variables['time_series_600'])
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
        self.spec_x.opener(data)
        self.spec_y.opener(data)
        for step in range(nsteps):
            self.step=step
            helper.update(data,step)
            self.stat_1d.opener2(data,step)
            self.samp_1d.opener2(data,step)
            self.cross_xz.opener2(data,step)
            self.cross_yz.opener2(data,step)
            self.cross_xy.opener2(data,step)
            self.stat_xz.opener2(data,step)
            self.stat_yz.opener2(data,step)
            self.stat_int.opener2(data,step)
            self.stat_dom.opener2(data,step)
            self.spec_x.opener2(data,step)
            self.spec_y.opener2(data,step)
            self.calc_masks()
            timestep=self.data.variables['time_series_600'][self.step]
            self.clouds=statgroupclouds("pdfdata/pdfdata."+sysconfig.exper+".%05d.%03d.nc"%(timestep,sysconfig.filenumber),'3d in-cloud variable fields at %05d seconds' %timestep)
            self.clouds.opener(data)
            self.clouds.opener2(data,step)
            self.processor()
            self.clouds.closer()
        self.stat_1d.closer()
        self.samp_1d.closer()
        self.cross_xz.closer()
        self.cross_yz.closer()
        self.cross_xy.closer()
        self.stat_xz.closer()
        self.stat_yz.closer()
        self.stat_int.closer()
        self.stat_dom.closer()
        self.spec_x.closer()
        self.spec_y.closer()
    def stat_var(self,var,field):
        self.stat_1d.put_make_var(var,mean_2d(field))      
        self.stat_xz.put_make_var(var,mean_xz(field))
        self.stat_yz.put_make_var(var,mean_yz(field))
        self.masked_var(var,field)
    def process_var(self,var,field):
        self.stat_1d.put_make_var(var,mean_2d(field))
        self.cross_xz.put_make_var(var,take(field,outputconfig.ysel,axis=1))
        self.cross_yz.put_make_var(var,take(field,outputconfig.xsel,axis=0)) 
        self.cross_xy.put_make_var(var,take(field,outputconfig.zsel,axis=2))           
        self.stat_xz.put_make_var(var,mean_xz(field))
        self.stat_yz.put_make_var(var,mean_yz(field))
        if 'makespectra' in allfields[var]:
            self.spec_x.put_make_var(var,field)
            self.spec_y.put_make_var(var,field)
        self.masked_var(var,field)
        #if var in progfields:
        if var in progfields or var in ['QT','THV']:
            dumpmasks=['tr_cc','tr_sh','br_cc','br_sh','bl_cc','bl_sh','tl_cc','tl_sh']
            if not outputconfig.lsamp:
                return
            elif 'xe' in allfields[var]:
                for mask in dumpmasks:
                    temp=1.0*field #force copy
                    whereinf=(self.masks[mask].fieldxe==0);
                    temp[whereinf] = nan
                    self.clouds.put_make_sampvar(var,temp,mask)          
            elif 'ye' in allfields[var]:
                for mask in dumpmasks:
                    temp=1.0*field #force copy
                    whereinf=(self.masks[mask].fieldye==0);
                    temp[whereinf] = nan
                    self.clouds.put_make_sampvar(var,temp,mask)          
            elif 'ze' in allfields[var]:
                for mask in dumpmasks:
                    temp=1.0*field #force copy
                    whereinf=(self.masks[mask].fieldze==0);
                    temp[whereinf] = nan
                    self.clouds.put_make_sampvar(var,temp,mask)          
            else:
                for mask in dumpmasks:
                    temp=1.0*field #force copy
                    whereinf=(self.masks[mask].field==0);
                    temp[whereinf] = nan
                    self.clouds.put_make_sampvar(var,temp,mask)  

    def ref_var(self,var,field):
        self.stat_1d.put_make_var(var,field)
    def masked_var(self,var,field):
        if not outputconfig.lsamp:
            return
        elif 'xe' in allfields[var]:
            for mask in self.masks.keys():         
                m1=field*self.masks[mask].fieldxe
                m2=mean_2d(m1)*self.masks[mask].ifracxe
                self.samp_1d.put_make_sampvar(var,m2,mask)           
        elif 'ye' in allfields[var]:
            for mask in self.masks.keys():         
                m1=field*self.masks[mask].fieldye
                m2=mean_2d(m1)*self.masks[mask].ifracye
                self.samp_1d.put_make_sampvar(var,m2,mask)  
        elif 'ze' in allfields[var]:
            for mask in self.masks.keys():         
                m1=field*self.masks[mask].fieldze
                m2=mean_2d(m1)*self.masks[mask].ifracze
                self.samp_1d.put_make_sampvar(var,m2,mask)                      
        else:
            for mask in self.masks.keys():      
                m1=field*self.masks[mask].field
                m2=mean_2d(m1)*self.masks[mask].ifrac
                self.samp_1d.put_make_sampvar(var,m2,mask)
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
    def integrate_rho_ze(self,field):
        return sum(0.25*(field[:,:,1:]+field[:,:,:-1])*((self.helper.ze[1:]-self.helper.ze[0:-1])*(self.helper.rhon[1:]+self.helper.rhon[0:-1]))[None,None,:],axis=2,dtype=numpy.float64)      
    def integrate_rho_zc(self,field):
        return sum(field[:,:,1:]*((self.helper.zc[1:]-self.helper.ze[0:-1])*self.helper.rhon[1:])[None,None,:],axis=2,dtype=numpy.float64)      
    def hvort(self,u,v):
        dvdx=(v[:,:,:]-vstack((v[-1,:,:][None,:,:],v[:-1,:,:])))/(self.helper.xe[1]-self.helper.xe[0])
        dudy=(u[:,:,:]-hstack((u[:,-1,:][:,None,:],u[:,:-1,:])))/(self.helper.ye[1]-self.helper.ye[0])
        return dvdx-dudy                            
    def vvort(self,u,v,w):
        dwdx=(w[:,:,:]-vstack((w[-1,:,:][None,:,:],w[:-1,:,:])))/(self.helper.xe[1]-self.helper.xe[0])
        dwdy=(w[:,:,:]-hstack((w[:,-1,:][:,None,:],w[:,:-1,:])))/(self.helper.ye[1]-self.helper.ye[0])
        dvdz=dstack(((v[:,:,1:]-v[:,:,:-1])/(self.helper.zc[1:]-self.helper.zc[0:-1])[None,None,:],0.0*v[:,:,-1]))
        dudz=dstack(((u[:,:,1:]-u[:,:,:-1])/(self.helper.zc[1:]-self.helper.zc[0:-1])[None,None,:],0.0*u[:,:,-1]))
        vvort=sqrt(ye_to_yc((dwdy-dvdz)**2)+xe_to_xc((dudz-dwdx)**2))
        vvort[:,:,0]=0.0                            
        return vvort
    def getooss(self,u,v,w):
        dvdx=(v[:,:,:]-vstack((v[-1,:,:][None,:,:],v[:-1,:,:])))/(self.helper.xe[1]-self.helper.xe[0])
        dudy=(u[:,:,:]-hstack((u[:,-1,:][:,None,:],u[:,:-1,:])))/(self.helper.ye[1]-self.helper.ye[0])
        oo1=0.5*(dvdx-dudy)
        ss1=0.5*(dvdx+dudy)
        del dvdx,dudy
        dwdy=(w[:,:,:]-hstack((w[:,-1,:][:,None,:],w[:,:-1,:])))/(self.helper.ye[1]-self.helper.ye[0])
        dvdz=dstack(((v[:,:,1:]-v[:,:,:-1])/(self.helper.zc[1:]-self.helper.zc[0:-1])[None,None,:],0.0*v[:,:,-1]))
        oo2=(dwdy-dvdz)
        ss2=(dwdy+dvdz)
        oo2[:,:,0]=0.0                            
        ss2[:,:,0]=0.0
        del dwdy,dvdz                             
        dwdx=(w[:,:,:]-vstack((w[-1,:,:][None,:,:],w[:-1,:,:])))/(self.helper.xe[1]-self.helper.xe[0])
        dudz=dstack(((u[:,:,1:]-u[:,:,:-1])/(self.helper.zc[1:]-self.helper.zc[0:-1])[None,None,:],0.0*u[:,:,-1]))
        oo3=0.5*(dudz-dwdx)
        ss3=0.5*(dudz+dwdx)
        oo3[:,:,0]=0.0                            
        ss3[:,:,0]=0.0 
        del dwdx,dudz 
        dudx=(u[:,:,:]-vstack((u[-1,:,:][None,:,:],u[:-1,:,:])))/(self.helper.xe[1]-self.helper.xe[0])
        dvdy=(v[:,:,:]-hstack((v[:,-1,:][:,None,:],v[:,:-1,:])))/(self.helper.ye[1]-self.helper.ye[0])
        dwdz=dstack((0.0*w[:,:,0],(w[:,:,1:]-w[:,:,:-1])/(self.helper.zc[1:]-self.helper.zc[0:-1])[None,None,:])) 
        oo=2.0*(ye_to_yc(xe_to_xc(oo1**2))+ze_to_zc(ye_to_yc(oo2**2))+ze_to_zc(xe_to_xc(oo3**2)))
        ss=2.0*(ye_to_yc(xe_to_xc(ss1**2))+ze_to_zc(ye_to_yc(ss2**2))+ze_to_zc(xe_to_xc(ss3**2))+dudx*dudx+dvdy*dvdy+dwdz*dwdz)                                          
        return oo,ss
