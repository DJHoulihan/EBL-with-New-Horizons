# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 12:02:37 2021

@author: Dennis Houlihan
"""

import numpy as np
from photutils import find_peaks
from glob import glob
from photutils import EPSFBuilder
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import NDData
from photutils.psf import extract_stars
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm


files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data = image[0].data


def stars(image_data,threshold):
    """ Find stars within an image.
    
    Other Functions Used
    --------------------
    find_peaks: astropy function that finds bright peaks in an image.
    sigma_clipped_stats: astropy function that calculates sigma-clipped 
                         statistics on the provided data.
    extract_stars: astropy function that extracts cutout images centered on 
                   stars defined in the input catalog(s).
        
    Parameters
    ----------
    image_data: an array of the image data.
    threshold: a data threshold value for the find_peaks function. It will 
               only detect sources above this threshold.
    
    Returns
    -------
    stars: an EPSFStars object of photutils.psf.epsf_stars module
    ** Make sure to go through and see if there are any misidentifications**
    """
    data = image_data
    
    peaks_tbl = find_peaks(data, threshold= 0.003)  
    peaks_tbl['peak_value'].info.format = '%.8g'  # for consistent table output  

    # This makes it so no peaks that are a certain distance from the edge are counted
    size = 15
    hsize = (size - 1) / 2
    x = peaks_tbl['x_peak']  
    y = peaks_tbl['y_peak']  
    mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))  

    # Making a final table of stars
    stars_tbl = Table()
    stars_tbl['x'] = x[mask]  
    stars_tbl['y'] = y[mask] 
    return stars_tbl



def PSF_maker(data,stars_tbl, oversample, maxiters):
    ''' Construct an ePSF from identified stars in an image.
    
    Other Functions Used
    --------------------
    find_peaks: astropy function that finds bright peaks in an image.
    sigma_clipped_stats: astropy function that calculates sigma-clipped 
                         statistics on the provided data.
    extract_stars: astropy function that extracts cutout images centered on 
                   stars defined in the input catalog(s).
    
    Parameters
    ----------
    data: the image data
    stars_tbl : a table of identified stars in the image
    oversample : The oversampling factor of the ePSF relative to the input 
                 stars along the x and y axes
    maxiters : The maximum number of iterations to perform.

    Returns
    -------
    epsf : The constructed ePSF.
    std : Standard deviation of the constructed ePSF.
    fitted_stars : The input stars with updated centers and fluxes derived from
                   fitting the output epsf.

    '''
      # Subtracting the background from the star cutouts used for the ePSF. 
    # In this case, I used the sigma-clipped median value as the background level
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)  
    data -= median_val 


    nddata = NDData(data=data)    # converting to NDData type - required for extract_stars function
    stars = extract_stars(nddata, stars_tbl, size=8)  # Extracting star cutouts from the table
   
    
    epsf_builder = EPSFBuilder(oversampling=oversample, maxiters=maxiters, progress_bar=True)  
    epsf, fitted_stars = epsf_builder(stars)
    std = np.std(epsf.data)
    
    # Plotting ePSF
    plt.figure()
    norm = simple_norm(epsf.data, 'log', percent=99.)
    plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
    plt.colorbar()
    
    return epsf, std, fitted_stars

data_new = data[40,:,:] + data[56,:,:] + data[71,:,:] + data[88,:,:] + data[103,:,:]
table = stars(data_new,0.003)
table.remove_rows([1,5])

PSF = PSF_maker(data_new,table,4,12)

plt.figure()
norm = simple_norm(PSF[0].data, 'log', percent=99.)
plt.imshow(PSF[0].data, norm=norm, origin='lower', cmap='viridis')
plt.colorbar()
