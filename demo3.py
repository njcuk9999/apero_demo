#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise 3: Calibrating a 2D image

- Step 1: Find the localization files (A and B)
- Step 2: Find the wavelength solution file
- Step 3: Find the preprocessed Proxima image
- Step 4: In python overplot the localization traces on the preprocessed image
- Step 5: In python add the wavelengths at 10 nm intervals to the step 4 plot
          (i.e. at …, 1200nm, 1210nm, 1220nm, 1230nm, …)

Created on 2023-02-08

@author: artigau, cook
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.interpolate import UnivariateSpline as IUSpline

# all relevant files
image_file = 'NIRPS_2023-01-20T08_42_08_941_pp.fits'
loc_A_file = 'NIRPS_2022-11-25T11_13_50_435_pp_loco_A.fits'
loc_B_file = 'NIRPS_2022-11-25T11_15_30_766_pp_loco_B.fits'
wave_A_file = 'NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits'

# read relevant data
loco_a = fits.getdata(loc_A_file)
loco_B = fits.getdata(loc_B_file)
wave_a = fits.getdata(wave_A_file)

# need to add a '4' as a constant. We trim the periphery of the pp image
#   calibrated images do not have the top/bottom/left/right pixels
image = fits.getdata(image_file)[4:-4, 4:-4][::-1, ::-1]

# find cuts for image display
vmin, vmax = np.nanpercentile(image, [5, 99])
# displaying image with cuts in greyscale
plt.imshow(image, vmin=vmin, vmax=vmax, aspect='auto', origin='lower',
           cmap="gray")
# define a pixel grid
xpix = np.arange(image.shape[0])

# wavelengths that will be plotted
vals = np.arange(950, 1950, 10)

# pixel index along the dispersion direction
index = np.arange(loco_a.shape[1])
for order_num in range(loco_a.shape[0]):
    # spline the wavelength solution onto the values that will be plotted
    #   to get the pixel positions of the wavelength values
    xvals = IUSpline(wave_a[order_num], xpix, ext=1)(vals)
    # get a mask of good y-pixel locations for this order (fiber A)
    #   (should not be outside the image
    good = (loco_a[order_num] > 0) & (loco_a[order_num] < loco_a.shape[1])
    plt.plot(index[good], loco_a[order_num][good], color='orange',
             linestyle='--', alpha=0.5)
    # only plot those pixel values not at zero
    good_vals = xvals != 0
    # don't plot if we have no values
    if np.sum(good_vals) != 0:
        # pixel positions that match good_vals
        xvals_plot = np.array(xvals[good_vals], dtype=int)
        # wavelength values that match good_vals
        vals_plot = vals[good_vals]
        # plot the localisation points as red dots
        plt.plot(xvals_plot, loco_a[order_num, xvals_plot], 'ro')
        # add text for each point (in yellow)
        for itxt in range(len(xvals_plot)):
            xplt = xvals_plot[itxt]
            yplt = loco_a[order_num, xvals_plot[itxt]]
            txt = '{:.0f}nm'.format(vals_plot[itxt])
            plt.text(xplt, yplt, txt, color='yellow')
    # get a mask of good y-pixel location for this order (fiber B)
    good = (loco_B[order_num] > 0) & (loco_B[order_num] < loco_B.shape[1])
    plt.plot(index[good], loco_B[order_num][good], color='green',
             linestyle='--', alpha=0.5)

plt.show()
