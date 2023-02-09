#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise 2: Looking at detector cosmetics

- Step 1: Find the DARK_FP RAMP file
- Step 2: In python remove horizontal striping
- Step 3: In python remove the “butterfly” pattern

Created on 2023-02-08

@author: artigau, cook
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# DARK+FP file
file = 'NIRPS_2023-02-06T12_54_00_118.fits'

# read the image
im = fits.getdata(file)

# find illuminated pixels
bright_pix = im > np.nanpercentile(im, 95)

# dummy image with pixels that are illuminated as NaNs
dark_pix = np.array(im)
dark_pix[bright_pix] = np.nan

# median structure along y axis
med = np.nanmedian(dark_pix, axis=1)

# plot the the horizontal striping
fig, ax = plt.subplots(nrows=1, ncols=2, sharey='all')
vmin, vmax = np.nanpercentile(im, [10, 85])
ax[0].imshow(im, vmin=vmin, vmax=vmax, aspect='auto')
ax[1].plot(med, np.arange(len(med)))
ax[1].set(ylabel='y pixel', xlabel='common noise', title='Median profile')
plt.show()

# replicate into 2d array and subtract
im -= np.repeat(med, im.shape[1]).reshape(im.shape)
# do the same on the dummy image with flux masked
dark_pix -= np.repeat(med, im.shape[1]).reshape(im.shape)
fits.writeto('correction1.fits', im, overwrite=True)
# save the image that has been masked
fits.writeto('correction1_dark_pix.fits', dark_pix, overwrite=True)

# cube that will contain individual amplifiers flipped along the readout
# direction
cube_amps = np.zeros([32, 4096, 128])
for iamp in range(32):
    amp = dark_pix[:, iamp * 128:iamp * 128 + 128]
    if (iamp % 2) == 1:
        # swap for the butterfly pattern
        amp = np.array(amp[:, ::-1])
    cube_amps[iamp] = amp

med_amp = np.nanmedian(cube_amps, axis=0)

# show the butterfly image
vmin, vmax = np.nanpercentile(med_amp, [1, 99])
plt.imshow(med_amp, aspect='auto', vmin=vmin, vmax=vmax)
plt.show()

# push the median amplifier values into the butterfly image
butterfly = np.zeros_like(im)
for iamp in range(32):
    if (iamp % 2) == 0:
        butterfly[:, iamp * 128:iamp * 128 + 128] = med_amp
    else:
        butterfly[:, iamp * 128:iamp * 128 + 128] = med_amp[:, ::-1]

# write to file to look at in ds9
fits.writeto('correction2.fits', im - butterfly, overwrite=True)
fits.writeto('butterfly.fits', butterfly, overwrite=True)
