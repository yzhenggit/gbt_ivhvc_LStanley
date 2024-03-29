import matplotlib.pyplot as plt
from yzGALFAHI.get_cubeinfo import get_cubeinfo
import numpy as np
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.wcs as wcs
import astropy.units as u
import astropy.constants as const
import sys


### first read in Lucas's data cube
ld_cube = 'data/mw_ta_cube_Lucas_Stanley.fits'
ld_data = fits.getdata(ld_cube)
ld_header = fits.getheader(ld_cube)
ld_data = ld_data[0]

# the header has 4D, and the fourth dimension is useless..., let's reduce that dimension
ld_header_3d = ld_header.copy()
ld_header_3d['NAXIS'] = 3
for ikey in ['NAXIS4', 'CTYPE4', 'CRVAL4', 'CRPIX4', 'CDELT4']:
    del ld_header_3d[ikey]

# sort out the ra, dec, and vlsr for Lucas's cube
ld_ra, ld_dec, ld_freq, ld_hdr_arrs = get_cubeinfo(ld_header_3d, returnHeader=True)
ld_hdr_2d = ld_hdr_arrs[0]
ld_wcs = wcs.WCS(ld_hdr_2d)

rest_freq = 1.420405751e9/u.s
rest_wave = const.c/rest_freq
ld_wave = const.c/(ld_freq/u.s)
ld_vlsr = ((ld_wave-rest_wave)/rest_wave*const.c).to(u.km/u.s).value

for pv in np.mgrid[-400:-99:10]:
    ### Now plot the GBT data 2D image:
    indv = np.argmin(abs(ld_vlsr-pv))

    figsize = (4.5, 6)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.1, 0.1, 0.75, 0.8], projection=ld_wcs)
    im = ax.imshow(ld_data[indv], origin='lower',
                   cmap=plt.cm.Blues, vmin=0.001, vmax=0.05)
    ax.coords[0].set_major_formatter('d.d')
    ax.coords[1].set_major_formatter('d.d')

    ra_max = 23.15
    ra_min = 21.8
    dec_min = 28.4
    dec_max = 31.1

    # set the boundaries of the plot
    x_min, y_min = ld_wcs.all_world2pix(ra_max, dec_min, 0)
    x_max, y_max = ld_wcs.all_world2pix(ra_min, dec_max, 0)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.set_title('vlsr=%.2f km/s\nMap saturated to 50mK'%(ld_vlsr[indv]), fontsize=14)
    ax.set_xlabel('RA (J2000)', fontsize=14)
    ax.set_ylabel('Dec (J2000)', fontsize=14)
    cb = fig.colorbar(im)
    cb.set_label('Tb (K)', fontsize=14)
    figfile = 'figs/LD_GBT_all_vchans/LD_GBT_vchan_%.2fkms.pdf'%(ld_vlsr[indv])
    fig.savefig(figfile)
    print(figfile)
    plt.close()
