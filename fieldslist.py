# STRUCTURE OF THIS FILE
# 'FIELDNAME':['description','units','arguments']
# arguments can currently relate to staggering, taking spectra

# prognostic variables at current time step
progfields={
'THETA':['potential temperature',u'K'],
'P':['pressure',u'Pa'],
'U':['wind speed in x direction',u'm s-1','xe','makespectra'],
'V':['wind speed in y direction',u'm s-1','ye','makespectra'],
'W':['wind speed in z direction',u'm s-1','ze','makespectra'],
'QV':['water vapor specific humidity',u'-'],
'QC':['cloud liquid water specific humidity',u'-'],
'QR':['rain specific humidity',u'-'],
'QS':['snow specific humidity',u'-'],
'QG':['graupel specific humidity',u'-'],
'QI':['cloud ice specific humidity',u'-'],
'NC':['cloud liquid number concentration',u'm-3'],
'NR':['rain number concentration',u'm-3'],
'NS':['snow number concentration',u'm-3'],
'NG':['graupel number concentration',u'm-3'],
'NI':['cloud ice number concentration',u'm-3'],
}

# prognostic variables at previous time step (necessary for leap-frog restart)
prevfields={
'THETAprev':['potential temperature, previous step',u'K',],
'Pprev':['pressure, previous step',u'Pa'],
'Uprev':['wind speed in x direction, previous step',u'm s-1','xe'],
'Vprev':['wind speed in y direction, previous step',u'm s-1','ye'],
'Wprev':['wind speed in z direction, previous step',u'm s-1','ze'],
'QVprev':['water vapor specific humidity, previous step',u'-'],
'QCprev':['cloud liquid water specific humidity, previous step',u'-'],
'QRprev':['rain specific humidity, previous step',u'-'],
'QSprev':['snow specific humidity, previous step',u'-'],
'QGprev':['graupel specific humidity, previous step',u'-'],
'QIprev':['cloud ice specific humidity, previous step',u'-'],
'NCprev':['cloud liquid number concentration, previous step',u'm-3'],
'NRprev':['rain number concentration, previous step',u'm-3'],
'NSprev':['snow number concentration, previous step',u'm-3'],
'NGprev':['graupel number concentration, previous step',u'm-3'],
'NIprev':['cloud ice number concentration, previous step',u'm-3'],
}

# derived variables (only at current time step)
derivedfields={
'T':['temperature',u'K'],
'RHO':['density (including hydrometeors, excluding pressure perturbation effects)',u'kg m-3'],
'THETAL':['linear liquid water potential temperature',u'K','makespectra'],
'THETARHO':['density potential temperature (incl hydrometeors)',u'K'],
'THETARHOX':['density potential temperature (excl hydrometeors)',u'K'],
'BUOY':['buoyancy (incl hydrometeors)',u'm s-2','spec'],
'BUOYX':['buoyancy (excl hydrometeors)',u'm s-2'],
'TMSE':['moist static energy temperature',u'K'],
'TLISE':['liquid ice static energy temperature (excl hydrometeors)',u'K'],
'QCI':['cloud water and ice specific humidity',u'-'],
'QP':['total precipitation specific humidity',u'-'],
'QT':['total water specific humidity (excl hydrometeors)',u'-','makespectra'],
'QSAT':['saturation specific humidity with respect to mixed phase',u'-'],
'QSATI':['saturation specific humidity with respect to ice phase',u'-'],
'RH':['relative humidity with respect to mixed phase',u'-'],
'RHI':['relative humidity with respect to ice phase',u'-'],
'UVAR':['variance of u',u'm2 s-2','xe'],
'VVAR':['variance of v',u'm2 s-2','ye'],
'WVAR':['variance of w',u'm2 s-2','ze'],
'THETAVAR':['variance of potential temperature',u'K2'],
'TMSEVAR':['variance of moist static energy temperature',u'K2'],
'TLISEVAR':['variance of liquid ice static energy temperature',u'K2'],
'THETALVAR':['variance of liquid water potential temperature',u'K2'],
'BUOYVAR':['variance of buoyancy',u'm2 s-4'],
'BUOYXVAR':['variance of buoyancy (exl hydrometeors)',u'm2 s-4'],
'QVVAR':['variance of water vapor',u'-'],
'QTVAR':['variance of total water',u'-'],
'RHOWBUOY':['anelastic buoyancy flux (incl hydrometeors)',u'kg m-1 s-3'],
'RHOWBUOYX':['anelastic buoyancy flux (excl hydrometeors)',u'kg m-1 s-3'],
'RHOWTHETA':['anelastic potential temperature flux',u'kg m-2 s-1 K'],
'RHOWQT':['anelastic total water flux (excl hydrometeors)',u'kg m-2 s-1'],
'RHOWMSE':['anelastic moist static energy temperature flux',u'kg m-2 s-1 K'],
'RHOWLISE':['anelastic liquid ice static energy temperature flux',u'kg m-2 s-1 K'],
'MASSFLX':['anelastic mass flux',u'kg m-2'],
'TKERES':['resolved TKE',u'm2 s-2'],
'TKESFS':['subfilter-scale TKE',u'm2 s-2'],
'KM':['eddy viscosity',u'm-1 s-1'],
'DP':['pressure perturbation with respect to slab mean',u'Pa'],
'DPDZ':['pressure gradient term',u'm s-2'],
'RHOWBMINP':['bouyancy minus pressure gradient term',u'kg m-1 s-3'],
'BMINP':['bouyancy minus pressure gradient term',u'm s-2'],
'HVORT':['horizontal component of the vorticity vector',u's-1'],
'VVORT':['vertical component of the vorticity vector',u's-1','ze'],
'BVWET2EXC':['moist BV frequency squared excess',u's-2','ze'],
'BG':['buoyancy generation',u'kg m-1 s-3','ze'],
}

# vertically integrated/maximum/minimum variables
intfields={
'VWP':['water wapor path',u'kg m-2'],
'CWP':['cloud liquid water path',u'kg m-2'],
'IWP':['cloud ice path',u'kg m-2'],
'TWP':['non-precipiating total water path',u'kg m-2'],
'RWP':['rain water path',u'kg m-2'],
'GWP':['graupel water path',u'kg m-2'],
'SWP':['snow path',u'kg m-2'],
'WMAX':['maximum vertical velocity in column',u'm s-1','max'],
'WMIN':['minimum vertical velocity in column',u'm s-1','min'],
'BUOYXMAX':['maximum buoyancy in column',u'm s-2','max'],
'BUOYXMIN':['minimum buoyancy in column',u'm s-2','min'],
'CLDTOP':['cloud top',u'm','max'],
'CLDW1TOP':['cloudy updraft (1 m/s) top',u'm','max'],
'MAXWIND':['maximum wind speed',u'm s-1','max'],
'MAXWINDHEIGHT':['height of maximum wind speed',u'm'],
'RHOWINT':['integrated vertical mass-flux',u'kg m2 s-2'],
'RHOUVINT':['integrated horizontal mass-flux',u'kg m2 s-2'],
}

fracfields={
'frac':['fraction',u'-'],
'fracxe':['fraction (at xe)',u'-','xe'],
'fracye':['fraction (at ye)',u'-','ye'],
'fracze':['fraction (at ze)',u'-','ze'],
}

# domain only variables
domfields={
'CC':['cloud cover',u'-'],
}

# variables with vertical dependence only
# note that cloud fractions etc are derived by introducing masks
vertfields={
'RHOREF':['anelastic reference density at p levels',u'kg m-3'],
'RHOREFH':['anelastic reference density at w levels',u'kg m-3','ze'],
'EXNREF':['anelastic exner function',u'-'],
'PREF':['anelastic reference pressure',u'-'],
'THETAREF':['anelastic reference temperature',u'-'],
'BVDRY2':['BV frequency squared',u's-2','ze']
}

