#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise 5: Telluric correction

- Step 1: Get the tcorr LBL file for Proxima
- Step 2: Get the associated recon file for the tcorr Proxima file.
- Step 3: Plot a histogram of the sigma away from the mean velocity
- Step 4: Plot the “trumpet plot” with the 3.5 sigma outliers and see if they match absorption features.

Created on 2023-02-08

@author: artigau, cook
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit


def gaussian(xvector, ew, amp):
    """
    Gaussian model
    :param xvector: x pixel positions
    :param ew: equivlent width (FWHM)
    :param amp: the amplitude of the gaussian

    :return: np.array, same size as xpix
    """
    # gaussian fit centered at zero
    yvector = np.exp(-0.5 * (xvector / ew) ** 2) * amp
    return yvector


def odd_ratio_mean(value, error, odd_ratio=1e-4, nmax=10,
                   conv_cut=1e-2):
    """
    Provide values and corresponding errors and compute a weighted mean
    advanced way to do the weighted mean and error on the weighted mean

    :param value: np.array (1D), value array
    :param error: np.array (1D), uncertainties for value array
    :param odd_ratio: float, the probability that the point is bad
                    Recommended value in Artigau et al. 2021 : f0 = 0.002
    :param nmax: int, maximum number of iterations to pass through
    :param conv_cut: float, the convergence cut criteria - how precise we have
                     to get

    :return: tuple, 1. the weighted mean, 2. the error on weighted mean
    """
    # deal with NaNs in value or error
    keep = np.isfinite(value) & np.isfinite(error)
    # deal with no finite values
    if np.sum(keep) == 0:
        return np.nan, np.nan
    # remove NaNs from arrays
    value, error = value[keep], error[keep]
    # work out some values to speed up loop
    error2 = error ** 2
    # placeholders for the "while" below
    guess_prev = np.inf
    # the 'guess' must be started as close as we possibly can to the actual
    # value. Starting beyond ~3 sigma (or whatever the odd_ratio implies)
    # would lead to the rejection of pretty much all points and would
    # completely mess the convergence of the loop
    guess = np.nanmedian(value)
    bulk_error = 1.0
    ite = 0
    # loop around until we do all required iterations
    while (np.abs(guess - guess_prev) / bulk_error > conv_cut) and (ite < nmax):
        # store the previous guess
        guess_prev = float(guess)
        # model points as gaussian weighted by likelihood of being a valid point
        # nearly but not exactly one for low-sigma values
        gfit = (1 - odd_ratio) * np.exp(-0.5 * ((value - guess) ** 2 / error2))
        # find the probability that a point is bad
        odd_bad = odd_ratio / (gfit + odd_ratio)
        # find the probability that a point is good
        odd_good = 1 - odd_bad
        # calculate the weights based on the probability of being good
        weights = odd_good / error2
        # update the guess based on the weights
        if np.sum(np.isfinite(weights)) == 0:
            guess = np.nan
        else:
            guess = np.nansum(value * weights) / np.nansum(weights)
            # work out the bulk error
            bulk_error = np.sqrt(1.0 / np.nansum(odd_good / error2))
        # keep track of the number of iterations
        ite += 1
    # return the guess and bulk error
    return guess, bulk_error


# relevant files
file_lblrv = 'NIRPS_2023-01-20T08_42_08_941_pp_e2dsff_tcorr_A_PROXIMA_PROXIMA_lbl.fits'
file_recon = 'NIRPS_2023-01-20T08_42_08_941_pp_s1d_v_recon_A.fits'

# reconstructed absorption
recon = Table.read(file_recon)

# read the LBL file
tbl = Table.read(file_lblrv)
# anything with <20 m/s per line is suspect. To be discussed with team
bad = tbl['sdv'] < 20
tbl['dv'][bad] = np.nan
tbl['sdv'][bad] = np.nan

# keep only valid entries in the table
tbl = tbl[np.isfinite(tbl['dv'])]

# get the mean velocity

# the best way is with a weighted mean
mu, sig = odd_ratio_mean(tbl['dv'], tbl['sdv'])
# but as a simple example we could just take a median and sig of good lines
"""
good_lines = tbl['sdv'] > 20 & tbl['sdv'] < 200
mu = np.nanmedian(tbl['dv'][good_lines])
sig = np.nanstd(tbl['dv'][good_lines])/np.sqrt(np.sum(good_lines) - 1)
"""

# plot the histogram of RV accuracy per line
v = np.array(tbl['sdv'])
plt.hist(v[np.isfinite(v)], range=[0, 1000], bins=1000)
plt.xlabel('RV accuracy per line [m/s]')
plt.ylabel('Number of lines')
plt.show()

# N sigma of each line relative to the mean velocity
nsig = np.array((tbl['dv'] - mu) / tbl['sdv'])
# distribution of N sigma
xh = plt.hist(nsig[np.isfinite(nsig)], bins=100)
yy = xh[0]
xx = (xh[1][:-1] + xh[1][1:]) / 2.0

# starting point for the curve fit
guess = [1.0, np.max(yy)]
# gaussian curve fit to the distribution
fit, _ = curve_fit(gaussian, xx, yy, guess=p0)
print('Dispersion relative to expected : {:.4f}'.format(fit[0]))
plt.plot(xx, gaussian(xx, *fit), 'r--')
plt.show()

# plot  the position of spurious points compared to telluric absorption
spurious = np.abs(nsig) > 3.5
fig, ax = plt.subplots(nrows=2, ncols=1, sharex='all')
# for each line, get the central position
ww = (tbl['WAVE_START'] + tbl['WAVE_END']) / 2.0
# get the width of each line
err_wave = (tbl['WAVE_END'] - tbl['WAVE_START']) / 2.0
ax[0].errorbar(ww, tbl['dv'], yerr=tbl['sdv'], xerr=err_wave, fmt='b.',
               alpha=0.1)
ax[0].errorbar(ww[spurious], tbl['dv'][spurious], yerr=tbl['sdv'][spurious],
               xerr=err_wave[spurious], fmt='r.', alpha=0.9)
# overplot the 1-sigma velocity error
ax[0].plot([np.min(tbl['WAVE_START']), np.max(tbl['WAVE_START'])],
           [mu, mu], 'r-')
ax[0].plot([np.min(tbl['WAVE_START']), np.max(tbl['WAVE_START'])],
           [mu - sig, mu - sig], 'r--')
ax[0].plot([np.min(tbl['WAVE_START']), np.max(tbl['WAVE_START'])],
           [mu + sig, mu + sig], 'r--')
ax[0].set(ylabel='Velocity')
# plot telluric absorption
ax[1].plot(recon['wavelength'], recon['flux'])
ax[1].set(xlabel='Wavelength [nm]', ylabel='Atmospheric transmission')
plt.show()
