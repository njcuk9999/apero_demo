# APERO Demos

The aim of these demos is to walk you through some of the key steps in the
data reduction using APERO.

- [Slides](https://docs.google.com/presentation/d/1fIolyXW1aPX5437BeNRkrmyukoeHmqudNnYj7YaGB04/edit?usp=sharing)
- [Example solutions](https://github.com/njcuk9999/apero_demo)

## Prerequisites for exercises 

- Python 3.9  (e.g. via [conda](https://docs.conda.io/en/latest/miniconda.html))
    ```
    conda create --name=apero-demo python=3.9
    conda activate apero-demo
    pip install -r requirements.txt
    ```


- DS9 ([download](https://sites.google.com/cfa.harvard.edu/saoimageds9)) 


- dfits and fitsort
  - python implementation: [download](https://astrom-tom.github.io/dfitspy/build/html/installation.html)
  - C implmentation: [download](https://github.com/granttremblay/eso_fits_tools)


- Download the file bundle: [download](https://www.astro.umontreal.ca/~artigau/apero_demo/apero_nirps_demo.tar)


## Exercise 1: Cube to RAMP: Correlated double sampling (CDS) 

- Step 1: Find the ramp and the cube for Proxima (HE)   HIERARCH ESO DPR TYPE = OBJECT,SKY  using dfits and fitsort (or python)
- Step 2: Load the cube in ds9
- Step 3: In DS9 play with the cube scaling (linear, log, histogram, min max, zscale etc)
- Step 4: In DS9 “Animate” the cube to see photons accumulating
- Step 5: In python create CDS “last minus first frame” from the cube. Express the resulting image in ADU/s
- Step 6: In DS9 compare to provided ramp image for Proxima

Example answer code: `demo1.py`

## Exercise 2: Looking at detector cosmetics

- Step 1: Find the DARK_FP RAMP file
- Step 2: In python remove horizontal striping 
- Step 3: In python remove the “butterfly” pattern

Example answer code: `demo2.py`

## Exercise 3: Calibrating a 2D image 

- Step 1: Find the localization files (A and B)
- Step 2: Find the wavelength solution file
- Step 3: Find the preprocessed Proxima image
- Step 4: In python overplot the localization traces on the preprocessed image
- Step 5: In python add the wavelengths at 10 nm intervals to the step 4 plot
          (i.e. at …, 1200nm, 1210nm, 1220nm, 1230nm, …)

Example answer code: `demo3.py`

## Exercise 4: Extraction, making the S1D

- Step 1: Find the e2dsff A file for Proxima
- Step 2: Find the associated blaze file using the e2dsff header
- Step 3: Get the wavelength solution from the e2dsff (header or ORDER_TABLE extension)
- Step 4: In python construct a 1D wavelength grid
- Step 5: In python calculate an S1D using the blaze and the flux

Example answer code: `demo4.py`

## Exercise 5: Telluric correction

- Step 1: Find the extracted s1d for Proxima (fiber A)
- Step 2: Find the telluric corrected s1d for Proxima (fiber A)
- Step 3: Find the associated recon and s1d template
- Step 4: Shift the s1d template using the BERV
- Step 5: Plot the flux (corrected/uncorrected) and template
- Step 6: Plot the normalized ratio of the corrected flux to template
- Step 7: Plot the transmission
- Bonus: Plot the OH line model from the pclean file (use fitsinfo to look at the extension names)

Example answer code: `demo5.py`

## Exercise 6: Radial velocities using the line-by-line approach

- Step 1: Get the tcorr LBL file for Proxima
- Step 2: Get the associated recon file for the tcorr Proxima file.
- Step 3: Plot a histogram of the sigma away from the mean velocity
- Step 4: Plot the “trumpet plot” with the 3.5 sigma outliers and see if they match absorption features.

Example answer code: `demo6.py`