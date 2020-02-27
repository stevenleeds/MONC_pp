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
        thetaref=self.gref('thref')
        delp=self.gv('p')
        buoyp=self.gv('buoy_p')
        buoyp=delp-buoyp
        qv=self.gv('q_vapour')
        qc=self.gv('q_cloud_liquid_mass')
        qtracer=self.gv('q_qfield_4')
        qi=0.0*qc
        qci=qc+qi
        pref=self.gref('prefn')
        p=delp+pref[None,None,:]
        exn=(pref/psfr)**(rd/cp)
        theta=thetaref[None,None,:]+deltheta
        del deltheta,delp
        t=theta*exn[None,None,:]
        self.process_var('U',u)      
        self.process_var('V',v)      
        self.process_var('W',w)
        self.process_var('QC',qc)
        self.process_var('QV',qv)
        qt=qci+qv
        self.process_var('QT',qt)
        self.process_var('TRACER',qtracer)
        self.process_var('THETA',theta)  
        self.process_var('P',p)
        self.process_var('BUOYP',buoyp)
        self.process_var('T',t)
        self.process_var('TMSE',t+(grav/cp)*self.helper.zc[None,None,:]+(rlvap/cp)*qv-((rlsub-rlvap)/cp)*qi)
        tmse=t+(grav/cp)*self.helper.zc[None,None,:]+(rlvap/cp)*qv-((rlsub-rlvap)/cp)*qi
        tlise=t+(grav/cp)*self.helper.zc[None,None,:]-(rlvap/cp)*qc-(rlsub/cp)*qi
        self.process_var('TLISE',tlise)
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
        self.ref_var('RHOREF',rhon)
        self.ref_var('THETAREF',thetaref)
        self.ref_var('RHOREFH',rho)
        thetarhox=theta*(1+(rvord-1)*qv-qci)
        del qci 
        self.process_var('THETARHOX',thetarhox)
        buoyx=grav*deviation_2d(thetarhox)/thetaref[None,None,:]
        self.process_var('BUOYX',buoyx)  
        self.stat_var('BUOYXVAR',buoyx*buoyx)      
        self.process_var('RHOWBUOYX',buoyx*w*rho[None,None,:])
        bvdry2=mean_2d((grav/thetaref[None,None,:])*dz(concatenate((thetarhox[:,:,1][:,:,None],thetarhox[:,:,1:]),axis=2),self.helper.zc))
        self.ref_var('BVDRY2',bvdry2)
        self.process_var('BVWET2EXC',bvwet2-bvdry2)
        self.process_var('BG',rho[None,None,:]*w*(bvwet2-bvdry2))
        del thetarhox,bvdry2,bvwet2
        rqtracer=rhon[None,None,:]*qtracer
        del qtracer
        qzeta=self.gv('q_qfield_3')
        self.process_var('QZETA',qzeta)
        wfull=ze_to_zc(w)
        updraft=(wfull*(wfull>0.))
        self.dom_var('TRACER_C',mean(rqtracer))
        self.dom_var('TRACER_CC',mean(rqtracer*rqtracer))
        self.dom_var('TRACER_CCC',mean(rqtracer*rqtracer*rqtracer))
        self.dom_var('TRACER_CB',mean(rqtracer*buoyx))
        self.dom_var('TRACER_CCB',mean(rqtracer*rqtracer*buoyx))
        self.dom_var('TRACER_CW',mean(rqtracer*wfull))
        self.dom_var('TRACER_CCW',mean(rqtracer*rqtracer*wfull))
        self.dom_var('TRACER_CZ',mean(rqtracer*self.helper.zc[None,None,:]))
        self.dom_var('TRACER_CCZ',mean(rqtracer*rqtracer*self.helper.zc[None,None,:]))
        self.dom_var('TRACER_CQ',mean(rqtracer*qt))
        self.dom_var('TRACER_CCQ',mean(rqtracer*rqtracer*qt))
        self.dom_var('TRACER_CH',mean(rqtracer*tmse))
        self.dom_var('TRACER_CCH',mean(rqtracer*rqtracer*tmse))
        self.dom_var('TRACER_CZETA',mean(rqtracer*qzeta))
        self.dom_var('TRACER_CCZETA',mean(rqtracer*rqtracer*qzeta))
        self.dom_var('TRACER_UC',mean(updraft*rqtracer))
        self.dom_var('TRACER_UCC',mean(updraft*rqtracer*rqtracer))
        self.dom_var('TRACER_UCCC',mean(updraft*rqtracer*rqtracer*rqtracer))
        self.dom_var('TRACER_UCB',mean(updraft*rqtracer*buoyx))
        self.dom_var('TRACER_UCCB',mean(updraft*rqtracer*rqtracer*buoyx))
        self.dom_var('TRACER_UCW',mean(updraft*rqtracer*wfull))
        self.dom_var('TRACER_UCCW',mean(updraft*rqtracer*rqtracer*wfull))
        self.dom_var('TRACER_UCZ',mean(updraft*rqtracer*self.helper.zc[None,None,:]))
        self.dom_var('TRACER_UCCZ',mean(updraft*rqtracer*rqtracer*self.helper.zc[None,None,:]))
        self.dom_var('TRACER_UCQ',mean(updraft*rqtracer*qt))
        self.dom_var('TRACER_UCCQ',mean(updraft*rqtracer*rqtracer*qt))
        self.dom_var('TRACER_UCH',mean(updraft*rqtracer*tmse))
        self.dom_var('TRACER_UCCH',mean(updraft*rqtracer*rqtracer*tmse))
        self.dom_var('TRACER_UCZETA',mean(updraft*rqtracer*qzeta))
        self.dom_var('TRACER_UCCZETA',mean(updraft*rqtracer*rqtracer*qzeta))
        self.dom_var('TRACER_U',mean(updraft))
        self.dom_var('TRACER_UB',mean(updraft*buoyx))
        self.dom_var('TRACER_UW',mean(updraft*wfull))
        self.dom_var('TRACER_UZ',mean(updraft*self.helper.zc[None,None,:]))
        self.dom_var('TRACER_UQ',mean(updraft*qt))
        self.dom_var('TRACER_UH',mean(updraft*tmse))
        self.dom_var('TRACER_UZETA',mean_2d(updraft*qzeta))
        self.ref_var('TRACERPROF_C',mean_2d(rqtracer))
        self.ref_var('TRACERPROF_CC',mean_2d(rqtracer*rqtracer))
        self.ref_var('TRACERPROF_CCC',mean_2d(rqtracer*rqtracer*rqtracer))
        self.ref_var('TRACERPROF_CB',mean_2d(rqtracer*buoyx))
        self.ref_var('TRACERPROF_CCB',mean_2d(rqtracer*rqtracer*buoyx))
        self.ref_var('TRACERPROF_CW',mean_2d(rqtracer*wfull))
        self.ref_var('TRACERPROF_CCW',mean_2d(rqtracer*rqtracer*wfull))
        self.ref_var('TRACERPROF_CZ',mean_2d(rqtracer*self.helper.zc[None,None,:]))
        self.ref_var('TRACERPROF_CCZ',mean_2d(rqtracer*rqtracer*self.helper.zc[None,None,:]))
        self.ref_var('TRACERPROF_CQ',mean_2d(rqtracer*qt))
        self.ref_var('TRACERPROF_CCQ',mean_2d(rqtracer*rqtracer*qt))
        self.ref_var('TRACERPROF_CH',mean_2d(rqtracer*tmse))
        self.ref_var('TRACERPROF_CCH',mean_2d(rqtracer*rqtracer*tmse))
        self.ref_var('TRACERPROF_CZETA',mean_2d(rqtracer*qzeta))
        self.ref_var('TRACERPROF_CCZETA',mean_2d(rqtracer*rqtracer*qzeta))
        self.ref_var('TRACERPROF_UC',mean_2d(updraft*rqtracer))
        self.ref_var('TRACERPROF_UCC',mean_2d(updraft*rqtracer*rqtracer))
        self.ref_var('TRACERPROF_UCCC',mean_2d(updraft*rqtracer*rqtracer*rqtracer))
        self.ref_var('TRACERPROF_UCB',mean_2d(updraft*rqtracer*buoyx))
        self.ref_var('TRACERPROF_UCCB',mean_2d(updraft*rqtracer*rqtracer*buoyx))
        self.ref_var('TRACERPROF_UCW',mean_2d(updraft*rqtracer*wfull))
        self.ref_var('TRACERPROF_UCCW',mean_2d(updraft*rqtracer*rqtracer*wfull))
        self.ref_var('TRACERPROF_UCZ',mean_2d(updraft*rqtracer*self.helper.zc[None,None,:]))
        self.ref_var('TRACERPROF_UCCZ',mean_2d(updraft*rqtracer*rqtracer*self.helper.zc[None,None,:]))
        self.ref_var('TRACERPROF_UCQ',mean_2d(updraft*rqtracer*qt))
        self.ref_var('TRACERPROF_UCCQ',mean_2d(updraft*rqtracer*rqtracer*qt))
        self.ref_var('TRACERPROF_UCH',mean_2d(updraft*rqtracer*tmse))
        self.ref_var('TRACERPROF_UCCH',mean_2d(updraft*rqtracer*rqtracer*tmse))
        self.ref_var('TRACERPROF_UCZETA',mean_2d(updraft*updraft*rqtracer*qzeta))
        self.ref_var('TRACERPROF_UCCZETA',mean_2d(updraft*rqtracer*rqtracer*qzeta))
        self.ref_var('TRACERPROF_U',mean_2d(updraft))
        self.ref_var('TRACERPROF_UB',mean_2d(updraft*buoyx))
        self.ref_var('TRACERPROF_UW',mean_2d(updraft*wfull))
        self.ref_var('TRACERPROF_UZ',mean_2d(updraft*self.helper.zc[None,None,:]))
        self.ref_var('TRACERPROF_UQ',mean_2d(updraft*qt))
        self.ref_var('TRACERPROF_UH',mean_2d(updraft*tmse))
        self.ref_var('TRACERPROF_UZETA',mean_2d(updraft*qzeta))
        del qzeta,wfull
        dp=deviation_2d(p)
        del p
        self.process_var('DP',dp)
        dpdz=(1./rho[None,None,:])*dz(dp,self.helper.zc)
        self.process_var('DPDZ',dpdz)
        dbuoyp=deviation_2d(buoyp)
        self.process_var('DBUOYP',dbuoyp)
        dbuoypdz=(1./rho[None,None,:])*dz(dbuoyp,self.helper.zc)
        self.process_var('DBUOYPDZ',dbuoypdz)
        self.dom_var('TRACER_CDPDZ',mean(rqtracer*dpdz))
        self.dom_var('TRACER_CCDPDZ',mean(rqtracer*rqtracer*dpdz))
        self.dom_var('TRACER_CDBUOYPDZ',mean(rqtracer*dbuoypdz))
        self.dom_var('TRACER_CCDBUOYPDZ',mean(rqtracer*rqtracer*dbuoypdz))
        self.dom_var('TRACER_UCDPDZ',mean(updraft*rqtracer*dpdz))
        self.dom_var('TRACER_UCCDPDZ',mean(updraft*rqtracer*rqtracer*dpdz))
        self.dom_var('TRACER_UCDBUOYPDZ',mean(updraft*rqtracer*dbuoypdz))
        self.dom_var('TRACER_UCCDBUOYPDZ',mean(updraft*rqtracer*rqtracer*dbuoypdz))
        self.ref_var('TRACERPROF_UCDPDZ',mean_2d(updraft*rqtracer*dpdz))
        self.ref_var('TRACERPROF_UCCDPDZ',mean_2d(updraft*rqtracer*rqtracer*dpdz))
        self.ref_var('TRACERPROF_UCDBUOYPDZ',mean_2d(updraft*rqtracer*dbuoypdz))
        self.ref_var('TRACERPROF_UCCDBUOYPDZ',mean_2d(updraft*rqtracer*rqtracer*dbuoypdz))
        self.ref_var('TRACERPROF_CDPDZ',mean_2d(rqtracer*dpdz))
        self.ref_var('TRACERPROF_CCDPDZ',mean_2d(rqtracer*rqtracer*dpdz))
        self.ref_var('TRACERPROF_CDBUOYPDZ',mean_2d(rqtracer*dbuoypdz))
        self.ref_var('TRACERPROF_CCDBUOYPDZ',mean_2d(rqtracer*rqtracer*dbuoypdz))
        del dp,dbuoyp,dbuoypdz,updraft
        self.process_var('BMINP',buoyx-dpdz)
        self.process_var('RHOWBMINP',(buoyx-dpdz)*w*rho[None,None,:])
        del dpdz,buoyx      
        thetal=theta-(rlvap/(cp*exn))*qc-(rlsub/(cp*exn))*qi
        self.process_var('THETAL',thetal) 
        del thetal
        qsat=qsa1/(0.01*pref[None,None,:]*exp(qsa2*(t-tk0c)/(t-qsa3))-qsa4) 
        qsati=qis1/(0.01*pref[None,None,:]*exp(qis2*(t-tk0c)/(t-qis3))-qis4)
        del t
        self.process_var('QSAT',qsat) 
        self.process_var('QSATI',qsati) 
        self.process_var('RH',qt/qsat)
        self.process_var('RHI',qt/qsati) 
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
        # additional scalars
        #for scalar_number in self.helper.svlist:
        #    self.process_var('SCALAR%03d'%scalar_number,self.gq('q',scalar_number)) 
        # height integrated variables 
        self.int_var('WMIN',self.helper.wmin)      
        self.int_var('WMAX',self.helper.wmax)
        self.int_var('CLDTOP',nanmax(self.helper.cld*self.helper.zc[None,None,:],axis=2))   
        self.int_var('CLDW1TOP',nanmax(self.helper.cld*(self.helper.wzc>1.0)*self.helper.zc[None,None,:],axis=2))
        self.int_var('VWP',self.integrate_rho_zc(qv))
        self.int_var('CWP',self.integrate_rho_zc(qc))
        self.int_var('IWP',self.integrate_rho_zc(qi))
        self.int_var('TWP',self.integrate_rho_zc(qt))
        del qv,qc,qi,qt
        # domain integrated variables
        vvort=self.vvort(u,v,w)
        self.process_var('VVORT',vvort)
        hvort=self.hvort(u,v)
        self.process_var('HVORT',hvort)
        tvort=sqrt(ye_to_yc(xe_to_xc(hvort**2))+ze_to_zc(vvort**2))
        self.process_var('TVORT',tvort)
        self.dom_var('TRACER_CTVORT',mean(rqtracer*tvort))
        self.dom_var('TRACER_CCTVORT',mean(rqtracer*rqtracer*tvort))
        self.ref_var('TRACERPROF_CTVORT',mean_2d(rqtracer*tvort))
        self.ref_var('TRACERPROF_CCTVORT',mean_2d(rqtracer*rqtracer*tvort))
        del hvort,vvort,tvort
        pp,qq,rr=self.vorticize(u,v,w)
        ccx,ccy,ccz=self.vorticize(pp,qq,rr)
        del pp,qq,rr
        self.process_var('CCX',ccx)
        self.process_var('CCY',ccy)
        self.process_var('CCZ',ccz)
        self.dom_var('TRACER_C_CCX',mean(rqtracer*ccx))
        self.dom_var('TRACER_CC_CCX',mean(rqtracer*rqtracer*ccx))
        self.ref_var('TRACERPROF_C_CCX',mean_2d(rqtracer*ccx))
        self.ref_var('TRACERPROF_CC_CCX',mean_2d(rqtracer*rqtracer*ccx))
        self.dom_var('TRACER_C_CCY',mean(rqtracer*ccy))
        self.dom_var('TRACER_CC_CCY',mean(rqtracer*rqtracer*ccy))
        self.ref_var('TRACERPROF_C_CCY',mean_2d(rqtracer*ccy))
        self.ref_var('TRACERPROF_CC_CCY',mean_2d(rqtracer*rqtracer*ccy))
        self.dom_var('TRACER_C_CCZ',mean(rqtracer*ccz))
        self.dom_var('TRACER_CC_CCZ',mean(rqtracer*rqtracer*ccz))
        self.ref_var('TRACERPROF_C_CCZ',mean_2d(rqtracer*ccz))
        self.ref_var('TRACERPROF_CC_CCZ',mean_2d(rqtracer*rqtracer*ccz))
        del ccx,ccy,ccz
        oo,ss=self.getooss(u,v,w)
        self.process_var('OO',oo)
        self.process_var('SS',ss)
        qq=0.5*(oo-ss)
        self.process_var('QCRITERION',qq)
        self.dom_var('TRACER_CQCRITERION',mean(rqtracer*qq))
        self.dom_var('TRACER_CCQCRITERION',mean(rqtracer*rqtracer*qq))
        self.ref_var('TRACERPROF_CQCRITERION',mean_2d(rqtracer*qq))
        self.ref_var('TRACERPROF_CCQCRITERION',mean_2d(rqtracer*rqtracer*qq))
        del oo,ss,qq,rqtracer
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
    make_onefile(sysconfig.filenumber-1)
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
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir+'/clouds')
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
    mkdir_p(multiplex_infrastructure.sysconfig.scratchdir+'/clouds')
    mkdir_p(multiplex_infrastructure.sysconfig.projectdir)
    runme()
    print('finished copying to project: '+case+' '+exper+' '+args.filenumber)
    writefinished('finished copying to project: '+args.case+' '+args.exper+' '+args.filenumber,multiplex_infrastructure.sysconfig.projectdir+'/finished.'+args.exper+'.'+args.filenumber+'.log')

