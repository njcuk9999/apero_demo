#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A crude slow all python version of dfits

Usage:

apero_dfits {filename}

apero_dfits '{glob command}' --fitsort {header keys}

i.e.

apero_dfits '*.fits' --fitsort OBJECT "HIERARCH ESO DPR TYPE"


Created on 2023-02-14

@author: artigau, cook
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import MinMaxInterval, ZScaleInterval
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1 import make_axes_locatable


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def get_args():
	"""
	Get arguments in a crude way (don't use argparse here)
	:return:
	"""
	filenames = sys.argv[1]
	return filenames


class Visualizer:
	def __init__(self, filename):
		# get data
		image = fits.getdata(filename)
		# set up figure
		self.fig = plt.figure()
		self.add_frames()
		# ---------------------------------------------------------------------
		self.current_interval = MinMaxInterval
		self.current_stretch = LinearStretch
		self.current_image = image
		self.current_file = filename

	def add_frames(self):
		self.main_frame = self.fig.add_subplot(1, 1, 1)
		self.fig.subplots_adjust(bottom=0.2)
		# ---------------------------------------------------------------------
		self.linear_stretch_axis = self.fig.add_axes([0.11, 0.05, 0.1, 0.075])
		self.linear_stretch_btn = Button(self.linear_stretch_axis, 'Linear')
		self.linear_stretch_btn.on_clicked(self.set_linear_strech)
		# ---------------------------------------------------------------------
		self.log_stretch_axis = self.fig.add_axes([0.22, 0.05, 0.1, 0.075])
		self.log_stretch_btn = Button(self.log_stretch_axis, 'Log')
		self.log_stretch_btn.on_clicked(self.set_log_strech)
		# ---------------------------------------------------------------------
		self.zscale_axis = self.fig.add_axes([0.33, 0.05, 0.1, 0.075])
		self.zscale_btn = Button(self.zscale_axis, 'ZScale')
		self.zscale_btn.on_clicked(self.set_z_scale)
		# ---------------------------------------------------------------------
		self.min_max_axis = self.fig.add_axes([0.44, 0.05, 0.1, 0.075])
		self.min_max_btn = Button(self.min_max_axis, 'MinMax')
		# ---------------------------------------------------------------------
		self.rotate_axis = self.fig.add_axes([0.55, 0.05, 0.1, 0.075])
		self.rotate_btn = Button(self.rotate_axis, 'Rotate90')
		self.rotate_btn.on_clicked(self.rotate_90)
		# ---------------------------------------------------------------------
		self.flipx_axis = self.fig.add_axes([0.66, 0.05, 0.1, 0.075])
		self.flipx_btn = Button(self.flipx_axis, 'FlipX')
		self.flipx_btn.on_clicked(self.flip_x)
		# ---------------------------------------------------------------------
		self.flipy_axis = self.fig.add_axes([0.77, 0.05, 0.1, 0.075])
		self.flipy_btn = Button(self.flipy_axis, 'FlipY')
		self.flipy_btn.on_clicked(self.flip_y)

	def imshow(self):
		# clear main frame
		self.fig.clf()
		# re-add frames
		self.add_frames()
		# get the new normalization
		norm = ImageNormalize(self.current_image,
		                      interval=self.current_interval(),
		                      stretch=self.current_stretch())
		# plot imshow
		im = self.main_frame.imshow(self.current_image, origin='lower',
		                            norm=norm, aspect='auto')
		divider = make_axes_locatable(self.main_frame)
		cax = divider.append_axes('right', size="5%", pad=0.1)
		self.fig.add_axes(cax)
		self.fig.colorbar(im, cax=cax, orientation="vertical")
		# set title
		plt.suptitle(self.current_file)
		# show figure
		plt.draw()

	def set_linear_strech(self, event):
		self.current_stretch = LinearStretch
		self.imshow()

	def set_log_strech(self, event):
		self.current_stretch = LogStretch
		self.imshow()

	def set_minmax_scale(self, event):
		self.current_interval = MinMaxInterval
		self.imshow()

	def set_z_scale(self, event):
		self.current_interval = ZScaleInterval
		self.imshow()

	def rotate_90(self, event):
		self.current_image = np.rot90(self.current_image)
		self.imshow()

	def flip_x(self, event):
		self.current_image = self.current_image[:, ::-1]
		self.imshow()

	def flip_y(self, event):
		self.current_image = self.current_image[::-1]
		self.imshow()

def main():
	"""
	Main function acts like dfits and fitsort (crude and slow)

	:return:
	"""
	# get filesnames
	filename = get_args()
	# make plot
	im = Visualizer(filename)
	im.imshow()
	plt.show()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
	# run main with no arguments (get from command line - sys.argv)
	ll = main()

# =============================================================================
# End of code
# =============================================================================
