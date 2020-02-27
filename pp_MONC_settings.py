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
                   
import os
import socket
import getpass
import ConfigParser
from json import loads
from numpy import array
## FILE PATHS

class opconfig():
    def __init__(self):
        pass
    def update(self,cfgfile):
        config = ConfigParser.RawConfigParser()
        config.read(cfgfile)
        self.lzlib=config.getboolean("outputs","lzlib")
        self.lcross_xy=config.getboolean("outputs","lcross_xy")
        self.lcross_xz=config.getboolean("outputs","lcross_xz")
        self.lcross_yz=config.getboolean("outputs","lcross_yz")
        self.ldiag_xz=config.getboolean("outputs","ldiag_xz")
        self.ldiag_yz=config.getboolean("outputs","ldiag_yz")
        self.ldiag_int=config.getboolean("outputs","ldiag_int")
        self.lsamp=config.getboolean("outputs","lsamplzlib")
        self.lspec=config.getboolean("outputs","lspec")
        self.lclouds=config.getboolean("outputs","lclouds")
        self.xsel=loads(config.get("domain","xsel"))
        self.ysel=loads(config.get("domain","ysel"))
        self.zsel=loads(config.get("domain","zsel"))
        self.dx=loads(config.get("domain","dx"))
        self.dy=loads(config.get("domain","dy"))
        self.spectralevelsbot=loads(config.get("domain","spectralevelsbot"))
        self.spectralevelstop=loads(config.get("domain","spectralevelstop"))
        self.nboundlines=config.getint("domain","nboundlines")
        self.xsel=array(self.xsel)-self.nboundlines
        self.ysel=array(self.ysel)-self.nboundlines

myusername=getpass.getuser()
hostname=socket.gethostname()
homedir = os.environ['HOME']

class sconfig():
    def __init__(self):
        pass
    def autoupdate(self,case,exper,overwrite):
        if 'see' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
            self.update('see.syscfg',case,exper,overwrite)
        elif 'arc2' in hostname and 'leeds' in hostname: ## UNIVERSITY OF LEEDS
            self.update('arc2.syscfg',case,exper,overwrite)
        elif 'xcm' in hostname: ## Met Office monsoon
            self.update('xcm.syscfg',case,exper,overwrite)
        elif 'nid' in hostname: ## Met Office monsoon
            self.update('xcm.syscfg',case,exper,overwrite)
        else: ## e.g. laptop
            self.update('laptop.syscfg',case,exper,overwrite)
    def do_replace(self,targ):        
        targ=targ.replace("$homedir",homedir)
        targ=targ.replace("$username",myusername)
        targ=targ.replace("$case",self.case)
        targ=targ.replace("$exper",self.exper)
        return targ
    def update(self,cfgfile,case,exper,overwrite):
        config = ConfigParser.RawConfigParser()
        config.read(cfgfile)
        self.case=case
        self.exper=exper
        self.overwrite=overwrite
        # fulldir: standard location of input directories
        # scratchdir: location of scratch space where initial postprocessing is done
        # localdir: location where final files are stored
        self.scratchdir=config.get("paths","scratchdir")
        self.projectdir=config.get("paths","projectdir")
        self.fulldir=config.get("paths","fulldir")
        self.scratchdir=self.do_replace(self.scratchdir)
        self.projectdir=self.do_replace(self.projectdir)
        self.fulldir=self.do_replace(self.fulldir)
