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

from multiplex_infrastructure import *
import multiplex_infrastructure
import argparse
   
numpy.seterr(invalid='ignore') # don't wine about nans

# class for postprocessing data
class dataprocessor(dataorganizer):
    # FUNCTIONS to refer to
    # process_var: process statistics and make cross sections
    # stat_var: process only statistics
    # int_var: this is a height integrated variable
    # ref_var: this is a reference variable
    def processor(self):
        u=self.gv('u_nogal')
        v=self.gv('v_nogal')
        w=self.gv('w')
        deltheta=self.gv('th')
        self.process_var('THPERT',deltheta)
        thetaref=self.gref('thref')
        delp=self.gv('p')
        buoyp=self.gv('buoy_p')
        buoyp=delp-buoyp
        qv=self.gv('q_vapour')
        qc=self.gv('q_cloud_liquid_mass')

        # Purity and radioactive tracers
        puritytracer=self.gv('q_purity_tracer')
        purity2=puritytracer*puritytracer
        radiotracer1=self.gv('q_radio1') #tau 30 mins
        radiotracer2=self.gv('q_radio2') #tau 5 mins

        # quadrant tracers
        trac_tl = self.gv('q_tracer_tl')
        trac_tr = self.gv('q_tracer_tr')
        trac_bl = self.gv('q_tracer_bl')
        trac_br = self.gv('q_tracer_br')

        qi=0.0*qc
        qci=qc+qi
        pref=self.gref('prefn')
        p=delp+pref[None,None,:]
        exn=(pref/psfr)**(rd/cp)
        theta=thetaref[None,None,:]+deltheta    
        del deltheta,delp
        t=theta*exn[None,None,:]
        self.process_var('T',t)    
        self.process_var('W',w)
        self.process_var('U',u)
        self.process_var('V',v)
        self.process_var('QC',qc)
        self.process_var('QV',qv)
        qt=qci+qv
        self.process_var('QT',qt)

        # Add in extra tracer fields
        self.process_var('PURITY',puritytracer)
        self.process_var('PURITY2',purity2)
        self.process_var('RADIO1',radiotracer1)
        self.process_var('RADIO2',radiotracer2)
        self.process_var('TRACER_TL',trac_tl)
        self.process_var('TRACER_TR',trac_tr)
        self.process_var('TRACER_BL',trac_bl)
        self.process_var('TRACER_BR',trac_br)

        #print 'shape of theta array is', np.shape(theta)
        self.process_var('THETA',theta)          
        tmse=t+(grav/cp)*self.helper.zc[None,None,:]+(rlvap/cp)*qv-((rlsub-rlvap)/cp)*qi
        tlise=t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*qc-(rlsub/cp)*qi
        zcplus=self.helper.zc+5.0
        prefplus=pref+5.0*dz_1d(pref,self.helper.zc)
        iexnplus=(prefplus/psfr)**-(rd/cp)
        (tvplus,qcplus)=fastmoistphysics(tlise,qt,zcplus,prefplus,force=self.force)
        del qcplus,zcplus,prefplus
        self.force=0
        zcminus=self.helper.zc-5.0
        prefminus=pref-5.0*dz_1d(pref,self.helper.zc)
        iexnminus=(prefminus/psfr)**-(rd/cp)
        (tvminus,qcminus)=fastmoistphysics(tlise,qt,zcminus,prefminus,force=self.force)
        del tlise
        bvwet2=(grav/thetaref[None,None,:])*(tvplus*iexnplus[None,None,:]-tvminus*iexnminus[None,None,:])/10.0
        del tvminus,tvplus,qcminus,zcminus,prefminus,iexnplus,iexnminus
        rhon=self.gref('rhon')       
        rho=self.gref('rho')
        self.ref_var('EXNREF',exn)
        self.ref_var('RHOREF',rhon) #at p levels
        self.ref_var('THETAREF',thetaref)
        ## ANNE ##
          #self.ref_var('RHOREFH',rho) #at w levels - anne
        # Make rho into a 3d array
        rho = np.array(rho)
        x_array = np.repeat(rho[np.newaxis,:], 160, axis=0)
        rho_array = np.repeat(x_array[:,np.newaxis,:], 160, axis=1)
        print 'shape of rho array is', np.shape(rho_array)
        self.process_var('RHOREFH', rho_array)
        #thetarhox=theta*(1+(rvord-1)*qv-qci) #################
        del qci 
        #rqtracer=rhon[None,None,:]*qtracer
        #del qtracer, p,    puritytracer, radiotracer1, radiotracer2
        del p,  puritytracer, purity2, radiotracer1, radiotracer2, trac_tl, trac_tr, trac_bl, trac_br

        thetal=theta-(rlvap/(cp*exn))*qc-(rlsub/(cp*exn))*qi
        self.process_var('THETAL', thetal)
        thetav=theta*(1+(0.61*qv)-qc) #virtual potential temperature
        del_thetav=deviation_2d(thetav)
        self.process_var('BUOY',del_thetav) #including some measure of buoyancy in the cross sections
        self.process_var('THV', thetav) #theta v

        del thetal,thetav
        qsat=qsa1/(0.01*pref[None,None,:]*exp(qsa2*(t-tk0c)/(t-qsa3))-qsa4) 
        qsati=qis1/(0.01*pref[None,None,:]*exp(qis2*(t-tk0c)/(t-qis3))-qis4)
        del t

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
        dqt=deviation_2d(qt)
        self.stat_var('QTVAR',dqt*dqt)
        self.stat_var('RHOWQT',dqt*w*rho[None,None,:])
        del dqt   
        
########################################################################
        # additional scalars
        for scalar_number in self.helper.svlist:
            self.process_var('SCALAR%03d'%scalar_number,self.gq('q',scalar_number)) 
        # height integrated variables 
        self.int_var('WMIN',self.helper.wmin)      
        self.int_var('WMAX',self.helper.wmax)
        self.int_var('CLDTOP',nanmax(self.helper.cld*self.helper.zc[None,None,:],axis=2))  
        self.int_var('CLDW1TOP',nanmax(self.helper.cld*(self.helper.wzc>1.0)*self.helper.zc[None,None,:],axis=2))
        self.int_var('VWP',self.integrate_rho_zc(qv))
        self.int_var('CWP',self.integrate_rho_zc(qc))
        self.int_var('IWP',self.integrate_rho_zc(qi))
        self.int_var('TWP',self.integrate_rho_zc(qt))


        #self.int_var('CLDTOP',nanmax(self.helper.cld*(-1e9)+self.helper.zc[None,None,:],axis=2))  
        self.int_var('CLDBASE',nanmin((1-self.helper.cld)*99999.9+self.helper.cld*self.helper.zc[None,None,:],axis=2))
        #self.int_var('CLDBASE',nanmin((1-self.helper.cld)*(-1e9)+self.helper.cld*self.helper.zc[None,None,:],axis=2))
        self.int_var('THV_MAXCLD',nanmax(self.helper.cld*del_thetav,axis=2)) #finding maximum perturbation of virtual potential temperature from the horizontal mean inside cloud

#############################################################################################
        del qv,qc,qi,qt,del_thetav
        
        #del rqtracer
        wspeed=sqrt(xe_to_xc(u*u)+ye_to_yc(v*v))
        self.int_var('MAXWIND',nanmax(wspeed,axis=2))
        self.int_var('MAXWINDHEIGHT',nanmax((wspeed>nanmax(wspeed,axis=2)[:,:,None]-1e-5)*self.helper.ze[None,None,:],axis=2))
        self.int_var('RHOUVINT',self.integrate_rho_zc(wspeed))
        self.int_var('RHOWINT',self.integrate_rho_ze(w))
        self.dom_var('CC',mean(nanmax(self.helper.cld,axis=2)))
        # produce areal coverage      
        if multiplex_infrastructure.outputconfig.lsamp:
            for mask in self.masks.keys():
               self.samp_1d.put_make_sampvar('frac',mean_2d(self.masks[mask].field),mask)                 
               self.samp_1d.put_make_sampvar('fracxe',mean_2d(self.masks[mask].fieldxe),mask)           
               self.samp_1d.put_make_sampvar('fracye',mean_2d(self.masks[mask].fieldye),mask)  
               self.samp_1d.put_make_sampvar('fracze',mean_2d(self.masks[mask].fieldze),mask)  

########### MAIN PROGRAM ########### 

def runme():
    update_variables()
    make_onefile(sysconfig.filenumber)
    current_processor=dataprocessor()
    process_onefile(current_processor)
    copy_files_to_project()

# ACTUALLY CALLS THE SCRIPT FROM THE COMMAND LINE
# Using the construction with
# if __name__ == "__main__"
# makes sure we can import the separate routines
# currently not used, but may be useful when analysing problems

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case")
    parser.add_argument("exper")
    parser.add_argument("filenumber")
    parser.add_argument ("-c", "--config", dest='config_file', default='default.cfg', type=str);
    parser.add_argument ("-s", "--sysconfig", dest='sysconfig_file', default=None, type=str);
    parser.add_argument ("-o", "--overwrite", dest='overwrite',action='store_true')
    args=parser.parse_args()
    multiplex_infrastructure.outputconfig.update(args.config_file)
    if args.sysconfig_file==None:
        multiplex_infrastructure.sysconfig.autoupdate(args.case,args.exper,args.filenumber,args.overwrite)
    else:
        multiplex_infrastructure.sysconfig.update(args.sysconfig_file,args.case,args.exper,args.filenumber,args.overwrite)
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir)
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir+'/pdfdata')
    mkdir_p(multiplex_infrastructure.sysconfig.projectdir)
    runme()
    print('finished copying to project: '+args.case+' '+args.exper+' '+args.filenumber)
    writefinished('finished copying to project: '+args.case+' '+args.exper+' '+args.filenumber,multiplex_infrastructure.sysconfig.projectdir+'/finished.'+args.exper+'.'+args.filenumber+'.log')

# another interface for profiling
def runprof(case,exper,filenumber,cfgfile='default.cfg',sysfile=None,overwrite=True):
    multiplex_infrastructure.outputconfig.update(cfgfile)
    if sysfile==None:
        multiplex_infrastructure.sysconfig.autoupdate(case,exper,filenumber,overwrite)
    else:
        multiplex_infrastructure.sysconfig.update(sysfile,case,exper,filenumber,overwrite)
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir)
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir+'/pdfdata')
    mkdir_p(multiplex_infrastructure.sysconfig.projectdir)
    runme()
    print('finished copying to project: '+case+' '+exper+' '+args.filenumber)
    writefinished('finished copying to project: '+args.case+' '+args.exper+' '+args.filenumber,multiplex_infrastructure.sysconfig.projectdir+'/finished.'+args.exper+'.'+args.filenumber+'.log')

