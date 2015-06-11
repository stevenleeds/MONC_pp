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

# TODO
# ADD FURTHER VARIABLES
# IN PARTICULAR HYDROMETEOR BUDGETS
# DISCUSS MOMENTUM, ENERGY, VORTICITY, HELICITY AND "PRESSURE LAPLACIAN" BUDGET TOOLS
# ADD AN OPTION PARSER?
# TAKE STAGGERING/WEIGHTING INTO ACCOUNT TO PRODUCE SPECTRA AT THE SAME LEVELS FOR W AND OTHER VARIABLES?
# FIND OUT IF SURFACE VALUES ARE CORRECTLY DEALT WITH

## IMPORTS

from pp_MONC_infrastructure import *
import pp_MONC_infrastructure
   
numpy.seterr(invalid='ignore') # don't wine about nans

# class for postprocessing data
class dataprocessor(dataorganizer):
    def processor(self):
        u=self.gv('u')
        v=self.gv('v')
        w=self.gv('w')
        deltheta=self.gv('th')
        thetaref=self.gref('thref')
        delp=self.gv('p')
        qv=self.gq('q',nqv)
        qc=self.gq('q',nqc)
        qi=self.gq('q',nqi)
	qci=qc+qi
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
        self.process_var('QT',qci+qv)
        self.process_var('THETA',theta)  
        self.process_var('P',p)
        self.process_var('T',t)
        self.process_var('TMSE',t+(grav/cp)*self.helper.zc[None,None,:]+(rlvap/cp)*qv-((rlsub-rlvap)/cp)*qi)
        self.process_var('TLISE',t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*qc-(rlsub/cp)*qi)
        rhon=self.gref('rhon')       
        rho=self.gref('rho')
        self.ref_var('EXNREF',exn)
        self.ref_var('RHOREF',rhon)
        self.ref_var('THETAREF',thetaref)
        self.ref_var('RHOREFH',rho)
        thetarhox=theta*(1+(rvord-1)*qv-qci)
        self.process_var('THETARHOX',thetarhox)
        buoyx=grav*deviation_2d(thetarhox)/thetaref[None,None,:]
        self.process_var('BUOYX',buoyx)  
        dbuoyx=deviation_2d(buoyx)
        self.stat_var('BUOYXVAR',dbuoyx*dbuoyx)      
        del dbuoyx
        del thetarhox
        dp=deviation_2d(p)
        self.process_var('DP',dp)
        zeroxy=0.0*dp[:,:,0][:,:,None]
        pgrad=(1./rho[None,None,:])*(concatenate((zeroxy,(dp[:,:,2:]-dp[:,:,:-2])/(self.helper.zc[None,None,2:]-self.helper.zc[None,None,:-2]),zeroxy),axis=2))
        self.process_var('PGRAD',pgrad)
        del dp
        self.process_var('BMINP',concatenate((zeroxy,buoyx[:,:,1:-1]-pgrad[:,:,1:-1],zeroxy),axis=2))
        del pgrad,buoyx      
        thetal=theta-(rlvap/(cp*exn))*qc-(rlsub/(cp*exn))*qi
        self.process_var('THETAL',thetal) 
        del thetal
        qsat=qsa1/(p*exp(qsa2*(t-tk0c)/(t-qsa3))-qsa4) 
        qsati=qis1/(p*exp(qis2*(t-tk0c)/(t-qis3))-qis4) 
        self.process_var('QSAT',qsat) 
        self.process_var('QSATI',qsati) 
        self.process_var('RH',(qci+qv)/qsat)
        self.process_var('RHI',(qci+qv)/qsati) 
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
        # height integrated variables 
        self.int_var('WMIN',self.helper.wmin)      
        self.int_var('WMAX',self.helper.wmax)
        self.int_var('CLDTOP',nanmax(self.helper.cld*self.helper.zc[None,None,:],axis=2))   
        self.int_var('CLDW1TOP',nanmax(self.helper.cld*(self.helper.wzc>1.0)*self.helper.zc[None,None,:],axis=2))
        # domain integrated variables 
        self.dom_var('CC',mean_2d(nanmax(self.helper.cld,axis=2)))      
        if lsamp:
            for mask in self.masks.keys():
               self.samp_1d.put_make_sampvar('frac',mean_2d(self.masks[mask].field),mask)                 
               self.samp_1d.put_make_sampvar('fracxe',mean_2d(self.masks[mask].fieldxe),mask)           
               self.samp_1d.put_make_sampvar('fracye',mean_2d(self.masks[mask].fieldye),mask)  
               self.samp_1d.put_make_sampvar('fracze',mean_2d(self.masks[mask].fieldze),mask)  

########### MAIN PROGRAM ########### 

def runme():
    update_variables()
    make_filelist()
    current_processor=dataprocessor()
    process_3doutput(current_processor)
    copy_files_to_project()

# currently not used, but may be useful when analysing problems
if loadmpl:
    import matplotlib
    matplotlib.use('agg') # first define agg output, then continue to load rest of matplotlib and pyplot 
    from matplotlib import *
    from matplotlib.pyplot import *

# ACTUALLY CALLS THE SCRIPT FROM THE COMMAND LINE
# Using the construction with
# if __name__ == "__main__"
# makes sure we can import the separate routines
 
if __name__ == "__main__":
    incase=sys.argv[1]
    inexper=sys.argv[2]
    pp_MONC_infrastructure.case=incase
    pp_MONC_infrastructure.exper=inexper
    pp_MONC_infrastructure.fulldir=fullbase+incase+'/output/'
    pp_MONC_infrastructure.scratchdir=scratchbase+incase+'/'
    pp_MONC_infrastructure.projectdir=projectbase+incase+'/'
    mkdir_p(pp_MONC_infrastructure.scratchdir)
    mkdir_p(pp_MONC_infrastructure.scratchdir+'/clouds')
    mkdir_p(pp_MONC_infrastructure.projectdir)
    runme()
