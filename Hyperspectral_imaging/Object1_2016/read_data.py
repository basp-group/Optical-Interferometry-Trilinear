
import pyfits
import scipy.io as sio

hdulist = pyfits.open('Object1-v3-2.oifits')

tbdata = hdulist[6].data

tbd = tbdata.field('T3AMP')
sio.savemat("t3amp.mat",{'t3amp':tbd})

tbd = tbdata.field('T3AMPERR')
sio.savemat("t3amperr.mat",{'t3amperr':tbd})

tbd = tbdata.field('T3PHI')
sio.savemat("t3phi.mat",{'t3phi':tbd})

tbd = tbdata.field('T3PHIERR')
sio.savemat("t3phierr.mat",{'t3phierr':tbd})

tbd = tbdata.field('U1COORD')
sio.savemat("u1coord.mat",{'u1coord':tbd})

tbd = tbdata.field('V1COORD')
sio.savemat("v1coord.mat",{'v1coord':tbd})

tbd = tbdata.field('U2COORD')
sio.savemat("u2coord.mat",{'u2coord':tbd})

tbd = tbdata.field('V2COORD')
sio.savemat("v2coord.mat",{'v2coord':tbd})


tbdata = hdulist[5].data
tbd = tbdata.field('VIS2DATA')
sio.savemat("vis2data.mat",{'vis2data':tbd})

tbd = tbdata.field('VIS2ERR')
sio.savemat("vis2err.mat",{'vis2err':tbd})

tbd = tbdata.field('UCOORD')
sio.savemat("ucoord.mat",{'ucoord':tbd})

tbd = tbdata.field('VCOORD')
sio.savemat("vcoord.mat",{'vcoord':tbd})


tbdata = hdulist[7].data
tbd = tbdata.field('FLUXDATA')
sio.savemat("fluxdata.mat",{'fluxdata':tbd})

tbd = tbdata.field('FLUXERR')
sio.savemat("fluxerr.mat",{'fluxerr':tbd})

tbdata = hdulist[3].data
tbd = tbdata.field('EFF_WAVE')
sio.savemat("eff_wave.mat",{'eff_wave':tbd})

tbd = tbdata.field('EFF_BAND')
sio.savemat("eff_band.mat",{'eff_band':tbd})

hdulist.close()
