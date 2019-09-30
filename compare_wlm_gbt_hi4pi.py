import matplotlib.pyplot as plt
from yzGALFAHI.get_cubeinfo import get_cubeinfo
import numpy as np
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.wcs as wcs
import astropy.units as u
import astropy.constants as const

## first read in the data of GBT
gbtdir = '/Users/Yong/Dropbox/GitRepo/Zheng19_WLMCGM/data/GBT'
gbt_files = ['GBT_PHL2525_1.txt', 'GBT_PHL2525_2.txt', 'GBT_PHL2525_3.txt',
             'GBT_PHL2525_4.txt', 'GBT_PHL2525_5.txt', 'GBT_PHL2525_6.txt',
             'GBT_PHL2525_7.txt', 'GBT_PHL2525_8.txt', 'GBT_PHL2525_9.txt',
             'GBT_PHL2525_10.txt', 'GBT_PHL2525_11.txt', 'GBT_PHL2525_12.txt',
              'GBT_PHL2525_13.txt']

# get coordinates and spectrum from each file
gbt_ra = np.zeros(len(gbt_files))
gbt_dec = np.zeros(len(gbt_files))
gbt_vlsr = []
gbt_flux = []
for i, ifile in enumerate(gbt_files):
    # coordinates
    f = open(gbtdir+'/'+ifile)
    for line in f:
        if line[0] == 'c': break
    f.close()

    aline = line.replace('\n', '')
    glon = np.float(aline.split(' ')[-2])
    glat = np.float(aline.split(' ')[-1])

    coord = SkyCoord(glon, glat, frame='galactic', unit=(u.deg, u.deg))
    gbt_ra[i] = coord.icrs.ra.deg
    gbt_dec[i] = coord.icrs.dec.deg

    # spectrum
    tb = Table.read(gbtdir+'/'+ifile, format='ascii', data_start=4)
    gbt_vlsr.append(tb['col1'])
    gbt_flux.append(tb['col2'])
    print(ifile, gbt_ra[i], gbt_dec[i], coord.galactic.l.degree, coord.galactic.b.degree)
gbt_ra[gbt_ra>10] = gbt_ra[gbt_ra>10]-360

## then read in the HI4PI data
# Get HI4PI cube:
hi4pi_cube = '/Users/Yong/Dropbox/databucket/HI4PI_cubes/CAR_D01.fits'
hi4pi_data = fits.getdata(hi4pi_cube)
hi4pi_hdr = fits.getheader(hi4pi_cube)

cube_ra, cube_dec, cube_vlsr, hdr_arrs = get_cubeinfo(hi4pi_hdr, returnHeader=True)
hi4pi_hdr2d = hdr_arrs[0]
hi4pi_wcs = wcs.WCS(hi4pi_hdr2d)

# now compare the flux of GBT and HI4PI
for i in range(gbt_ra.size):
    # i = 0
    ira = gbt_ra[i]
    idec = gbt_dec[i]
    try:
        px, py = hi4pi_wcs.all_world2pix(ira, idec, 0)
        px = int(px)
        py = int(py)
        hi4pi_ispec = hi4pi_data[:, py, px]
    except:
        print('yes')
        continue
    plt.figure()
    plt.plot(cube_vlsr, hi4pi_ispec, color='k', label='HI4PI')
    plt.plot(gbt_vlsr[i], gbt_flux[i], color='r', label='GBT/JL')
    plt.legend(fontsize=14)
    plt.xlim(-400, 400)
    # plt.ylim(-0.5, 6)
    plt.xlabel('Vlsr (km/s)', fontsize=14)
    plt.ylabel('Tb (K)', fontsize=14)
    plt.title('WLM_GBT_HI4PI_pointing%d.pdf'%(i+1))
    filename = 'figs/JL_GBT_WLM_HI4PI/WLM_GBT_HI4PI_pointing%d.pdf'%(i+1)
    plt.savefig(filename)

    plt.ylim(-0.1, 0.5)
    filename = 'figs/JL_GBT_WLM_HI4PI/WLM_GBT_HI4PI_pointing%d_zoom.pdf'%(i+1)
    plt.savefig(filename)
    print('Save to: ', filename)
    plt.close()
    # break
