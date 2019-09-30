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

# example:
# case 1: python compare_ivhvc_gbt_hi4pi.py
# case 2: python compare_ivhvc_gbt_hi4pi.py 23 31 -150 M33
#### first, decide which point to get spectra from
if len(sys.argv) == 5:
    pra = np.float(sys.argv[1]) # 22.8
    pdec = np.float(sys.argv[2]) # 30.9
    pv = np.float(sys.argv[3])
    point_tag = sys.argv[4]

elif len(sys.argv) > 1 and len(sys.argv) <5:
    print("Not enough arguments, try: ")
    print("python compare_ivhvc_gbt_hi4pi.py 23 31 -150 M33")
    sys.exit(0)
else:
    # if ra/dec are not given, then the default is M33's southern tip
    pra = 23
    pdec = 31
    pv = -220
    point_tag = 'M33'


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

### get spectra from the corresponding RA/DEC point
ld_px, ld_py = ld_wcs.all_world2pix(pra, pdec, 0)
ld_px = int(ld_px)
ld_py = int(ld_py)
ld_spec = ld_data[:, ld_py, ld_px]


### Now get HI4PI spectra
# Get HI4PI cube:
hi4pi_cube = '/Users/Yong/Dropbox/databucket/HI4PI_cubes/CAR_F02.fits'
hi4pi_data = fits.getdata(hi4pi_cube)
hi4pi_hdr = fits.getheader(hi4pi_cube)

hi4pi_ra, hi4pi_dec, hi4pi_vlsr, hi4pi_hdr_arrs = get_cubeinfo(hi4pi_hdr, returnHeader=True)
hi4pi_hdr2d = hi4pi_hdr_arrs[0]
hi4pi_wcs = wcs.WCS(hi4pi_hdr2d)

# get the flux
hi4pi_px, hi4pi_py = hi4pi_wcs.all_world2pix(pra, pdec, 0)
hi4pi_px = int(hi4pi_px)
hi4pi_py = int(hi4pi_py)

hi4pi_spec = hi4pi_data[:, hi4pi_py, hi4pi_px]

### Now plot the GBT data 2D image:
indv = np.argmin(abs(ld_vlsr-pv))

figsize = (4, 6)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes([0.1, 0.15, 0.8, 0.8], projection=ld_wcs)
ax.scatter([ld_px], [ld_py], color='r', s=40, facecolor='none', edgecolor='r', lw=2)
im = ax.imshow(ld_data[indv], origin='lower', cmap=plt.cm.Blues)
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

ax.set_title('vlsr=%.2f km/s'%(ld_vlsr[indv]), fontsize=14)
ax.set_xlabel('RA (J2000)', fontsize=14)
ax.set_ylabel('Dec (J2000)', fontsize=14)
cb = fig.colorbar(im)
cb.set_label('Tb (K)', fontsize=14)
fig.tight_layout()
fig.savefig('figs/%s_LD_GBT_RA%.2f_DEC%.2f_vchan.pdf'%(point_tag, pra, pdec, ld_vlsr[indv]))
plt.close()

## and plot the spectra
# now plot things
plt.figure()
plt.plot(hi4pi_vlsr, hi4pi_spec, color='k', label='HI4PI')
plt.plot(ld_vlsr,ld_spec, color='r', label='GBT/LD')
plt.title('GBT vs HI4PI, RA=%.2f, DEC=%.2f'%(pra, pdec))
plt.legend(fontsize=14)
plt.xlim(-500, 200)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Tb (km/s)')
plt.savefig('figs/%s_LD_GBT_RA%.2f_DEC%.2f_wHI4PI.pdf'%(point_tag, pra, pdec))

plt.ylim(-0.1, 0.5)
plt.savefig('figs/%s_LD_GBT_RA%.2f_DEC%.2f_wHI4PI_zoom.pdf'%(point_tag, pra, pdec))
plt.close()

print(pra, pdec, ld_px, ld_py)
