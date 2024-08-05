# EBL with New Horizons

This repository contains all the readings, code, and data I used in my Senior Capstone project at Rochester Institute of Technology. In this project I performed PSF photometry on the star Vega using images taken aboard NASA's New Horizons Spacecraft. 
The goal was to assess the spacecraft's ability to be used as to measure the Extragalactic Background Light (EBL) from an outer solar system vantage point. The abstract for my Capstone project is shown below. I have also attached my final Capstone paper as well as a document outlining the PSF photometry procedure.

## A Measurement of the Extragalactic Background Light

The extragalactic background light (EBL) is the summed emission from all sources outside of the
Milky Way Galaxy. An accurate measurement of the EBL can be used as a benchmark to test
whether there are any components in excess of the integrated light from galaxies. Any discrepancies
would imply the presence of sources of diuse extragalactic emission that are currently unaccounted
for, which is predicted to include light from the rst galaxies in the universe. We have studied
archival data from the Linear Etalon Imaging Spectral Array (LEISA) aboard the New Horizons
spacecraft to determine whether it could be used to measure the EBL at wavelengths between 1.2
to 2.5 $\mu$m. We have empirically determined LEISA's sensitivity to diffuse brightness to be $\delta\lambda$$I_\lambda$ =
2 x 10$^6$ nW/m/sr in 1.5 second exposures, which essentially precludes the possibility of reaching
EBL signals at the level of 10 nW/m/sr. These studies have allowed us to conclude that LEISA is
not even capable of generating useful upper limits on the amplitude of the EBL at near-infrared
wavelengths. We have also attempted to determine whether another New Horizons instrument, the
Multi-Spectral Visible Imaging Camera (MVIC), could produce a good measurement of the EBL.

## Programs Used
- Python (NumPy, SciPy, MatPlotLib, AstroPy, and other packages/libraries)
- SAOImage ds9: https://sites.google.com/cfa.harvard.edu/saoimageds9 
- Astrometry.net: https://nova.astrometry.net/upload

  ## Highlights:
<img src="https://github.com/user-attachments/assets/796075cb-1d61-44de-b658-a09b19299991" alt="2DGaussian" width="600" height="400">
<img src="https://github.com/user-attachments/assets/c549bc2e-152b-4906-b2a3-0da74f65b282" alt="ePSF" width="600" height="400">
<img src="https://github.com/user-attachments/assets/57a0eec3-eb84-4676-871f-ce27584c2149" alt="KLbestfit" width="600" height="400">
<img src="https://github.com/user-attachments/assets/ced79e90-4cf2-4597-8632-17544da02719" alt="residuals" width="600" height="400">

