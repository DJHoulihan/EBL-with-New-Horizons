# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 21:17:39 2021

@author: Dennis Houlihan
"""

import numpy as np
from photutils import find_peaks
from glob import glob
from astropy.visualization import ZScaleInterval 
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import NDData
from photutils.psf import extract_stars
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from photutils import EPSFBuilder
from scipy.interpolate import RegularGridInterpolator

files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Past Classes/RIT Undergrad/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data = image[0].data
# data = data[143,0:129,:]

# wdata = image[1].data
# interval = ZScaleInterval() 
# z = interval.get_limits(wdata[0])
# plt.imshow(wdata[0], origin='lower', cmap='gray')
# plt.colorbar()

def regrid(data, out_x, out_y):
    m = max(data.shape[0], data.shape[1])
    y = np.linspace(0, 1.0/m, data.shape[0])
    x = np.linspace(0, 1.0/m, data.shape[1])
    interpolating_function = RegularGridInterpolator((y, x), data)

    yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))

    return interpolating_function((xv, yv))

#%%

#                                         Generating an ePSF


# creating a summed image with different frames where vega is spread out
data_new = data[40,:,:] + data[56,:,:] + data[71,:,:] + data[88,:,:] + data[103,:,:]

interval = ZScaleInterval()    # these two lines are for plotting the image
z = interval.get_limits(data_new)  

#threshold = detect_threshold(data[185,:,:], 2)


# Finding peaks in the image and creating a table of them
peaks_tbl = find_peaks(data_new, threshold= 0.003)  
peaks_tbl['peak_value'].info.format = '%.8g'  # for consistent table output  
print(peaks_tbl)  


# This makes it so no peaks that are a certain distance from the edge are counted
# (ended up not really doing anything as nothing is close to the edge)
size = 15
hsize = (size - 1) / 2
x = peaks_tbl['x_peak']  
y = peaks_tbl['y_peak']  
mask = ((x > hsize) & (x < (data_new.shape[1] -1 - hsize)) & (y > hsize) & (y < (data_new.shape[0] -1 - hsize)))  


# Making a final table of stars
stars_tbl = Table()
stars_tbl['x'] = x[mask]  
stars_tbl['y'] = y[mask] 
print(stars_tbl)
stars_tbl.remove_rows([1,5])  # manually removing stars I saw were misidentified
print('New Table',stars_tbl)


    

# Subtracting the background from the star cutouts used for the ePSF. 
# In this case, I used the sigma-clipped median value as the background level
mean_val, median_val, std_val = sigma_clipped_stats(data_new, sigma=2.)  
data_new -= median_val 


nddata = NDData(data=data_new)    # converting to NDData type - required for extract_stars function
stars = extract_stars(nddata, stars_tbl, size=8)  # Extracting star cutouts from the table

# Scattering where the stars are on the image to make sure there are no misidentifications
plt.figure()
plt.imshow(data_new, cmap='gray',norm = None, vmin = z[0] , vmax = z[1])
plt.scatter(x[mask],y[mask], facecolors='none', edgecolors='r')

# Showing cutouts of each star with subplots
nrows = 2
ncols = 2
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

for i in range(nrows*ncols):
    norm = simple_norm(stars[i], 'log', percent=99.)
    ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
#plt.savefig('cutouts.pdf',bbox_inches = 'tight')
# Using epsf_builder function to create ePSF
epsf_builder = EPSFBuilder(oversampling=4, maxiters=12, progress_bar=True)  
epsf, fitted_stars = epsf_builder(stars)
std = np.std(epsf.data)

# Plotting ePSF
plt.figure()
norm = simple_norm(epsf.data, 'log', percent=99.)
plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
plt.colorbar()
#plt.savefig('ePSF.pdf',bbox_inches = 'tight')
#%%  Creating a window that will follow the star as it moves through the image

# I manually took the positions of Vega as it travelled through the image
# Each position is from a different frame 
y = np.array([62,66,65,63,61,58,61,64,67,63,60,58,61,64,66,64])
x = np.array([15,19,20,24,31,37,44,50,57,61,64,66,69,16,19,22])

# Calling the image 
files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Past Classes/RIT Undergrad/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data0 = image[0].data
data = data0[38:129,:,:] # cutting the data from frame 39 to frame 130

errdata0 = image[5].data
errdata = errdata0[38:129,:,:]
# Vega oscillates along the y-axis and travels up the x-axis. So I made the y length fixed
# by subtracting 20 pixels from its minimum y value and adding 20 pixels to its max y value.
# For the x-direction, I made the length 30 pixels centered around the first x value recorded. 
# Then for each frame, I add 1 pixel to the max and min x value of the stamp (called window) so it would
# follow the star as it moves up the x-axis.


wymin = min(y) - 20
wymax = max(y) + 20

wxmin = np.zeros(len(data),dtype = int)
wxmax = np.zeros(len(data), dtype = int)
middle = x[0] 
wxmin[0] = int(middle-15)
wxmax[0] = int(middle+15)


window = np.zeros(data.shape)
window = window[:,wxmin[0]:wxmax[0],wymin:wymax]
window[0] = data[0,wxmin[0]:wxmax[0],wymin:wymax]

errwindow = np.zeros(errdata.shape)
errwindow = errwindow[:,wxmin[0]:wxmax[0],wymin:wymax]
errwindow[0] = errdata[0,wxmin[0]:wxmax[0],wymin:wymax]


for i in range(len(data)-1):
   
    wxmin[i+1] = wxmin[i] + 1
    wxmax[i+1] = wxmax[i] + 1
    window[i+1] = data[i+1,wxmin[i+1]:wxmax[i+1], wymin:wymax]
    errwindow[i+1] = errdata[i+1,wxmin[i+1]:wxmax[i+1], wymin:wymax]

# i=55
# interval = ZScaleInterval() 
# z = interval.get_limits(window[i])
# plt.imshow(window[i], cmap='gray',norm = None, vmin = z[0] , vmax = z[1])

# Plotting each window to see if the stars are centered in each one
nrows = 6
ncols = 6
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

for i in range(nrows*ncols):
    norm = simple_norm(window[i], 'log', percent=99.)
    ax[i].imshow(window[i],norm=norm, origin='lower', cmap='viridis')
#plt.savefig('stamps.pdf',bbox_inches = 'tight')

#%%         Making Same window for the wavelength extension

# Repeating the same procedure for the wavlength extension of the image
y = np.array([62,66,65,63,61,58,61,64,67,63,60,58,61,64,66,64])
x = np.array([15,19,20,24,31,37,44,50,57,61,64,66,69,16,19,22])

wldata0 = image[1].data
wldata = wldata0[0,:,:] # convert from microns to Angstroms

middle = x[0]
wymin = min(y) - 20
wymax = max(y) + 20

wxmin = np.zeros(len(data),dtype = int)
wxmax = np.zeros(len(data), dtype = int)
wxmin[0] = int(middle-15)
wxmax[0] = int(middle+15)

wlwindow = np.zeros(window.shape)
wlwindow[0] = wldata[wxmin[0]:wxmax[0],wymin:wymax]

for i in range(len(wlwindow)-1):
   
    wxmin[i+1] = wxmin[i] + 1
    wxmax[i+1] = wxmax[i] + 1
    wlwindow[i+1,:,:] = wldata[wxmin[i+1]:wxmax[i+1], wymin:wymax]

nrows = 9
ncols = 9
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

for i in range(nrows*ncols):
    norm = simple_norm(wlwindow[i], 'log', percent=99.)
    ax[i].imshow(wlwindow[i],norm=norm, origin='lower', cmap='viridis')


#%%
nrows = 9
ncols = 9
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

for i in range(nrows*ncols):
    norm = simple_norm(scut[i,4:8,4:8], 'log', percent=99.)
    ax[i].imshow(scut[i,4:8,4:8],norm=norm, origin='lower', cmap='viridis')




#%% Using DAOStarFinder to get centroids of the star in each frame
from photutils import DAOStarFinder

sources = [] # creating empty list and will append sources into it
median = np.zeros((len(window)))
mean = np.zeros((len(window)))
std = np.zeros((len(window)))

# I fiddled with the inputs to DAOStarFinder to make sure only the star is detected,
# but there were still some misidentifications.
for i in range(len(window)):
    mean[i], median[i], std[i] = sigma_clipped_stats(window[i], sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=6.*std[i],roundlo = -2.0, roundhi = 2.0, sharplo = 0.2, sharphi = 1.1)
    sources.append(daofind(window[i]-median[i]))

# I manually went through the sources list and if there were more than one source detected
# in a list index, then I would go and delete the misidentifications.
sources[60].remove_row(0)
sources[25].remove_row(1)
sources[37].remove_row(1)
# sources[20].remove_row(1)
# sources[41].remove_row(1)
# sources[87].remove_row(0)
# sources[89].remove_rows([0,2])

# making an array of the stars positions as it's easier to work with (for me)
stars = np.zeros((len(sources),2))
for i in range(len(sources)):
    stars[i,0] = np.array(sources[i]['xcentroid'])
    stars[i,1] = np.array(sources[i]['ycentroid'])

#%% Creating Boxes around star centroids
positions = np.zeros((len(stars),2))
norm = np.zeros((len(stars),1))


nrows = 9
ncols = 9
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

xmin = np.zeros((len(stars)), dtype = int)
xmax = np.zeros((len(stars)), dtype = int)
ymin = np.zeros((len(stars)), dtype = int)
ymax = np.zeros((len(stars)), dtype = int)
#scut = np.zeros(window[0].shape)
radius = 4

scut = np.zeros((len(sources),radius*2,radius*2))
wcut = np.zeros((len(sources),radius*2,radius*2))
errcut = np.zeros((len(sources),radius*2,radius*2))
for i in range(len(sources)):
    positions[i] = np.transpose((stars[i,0], stars[i,1]))
 
    xmin[i] = int(round(positions[i,0]) - radius)
    xmax[i] = int(round(positions[i,0]) + radius)
    
    ymin[i] = int(round(positions[i,1]) - radius)
    ymax[i] = int(round(positions[i,1]) + radius)
    
    wcut[i] = wlwindow[i,ymin[i]:ymax[i], xmin[i]:xmax[i]]
    scut[i] = window[i,ymin[i]:ymax[i],xmin[i]:xmax[i]]
    errcut[i] = errwindow[i,ymin[i]:ymax[i],xmin[i]:xmax[i]]
for i in range(nrows*ncols):
    norm = simple_norm(window[i], 'log', percent=99.)
    ax[i].imshow(scut[i],norm=norm, origin='lower', cmap='viridis')
    #ax[i].scatter(positions[i,0],positions[i,1],facecolors='none', edgecolors='r')

plt.figure()
interval = ZScaleInterval() 
z = interval.get_limits(scut[0])
plt.imshow(scut[0], cmap='gray',norm = None, vmin = z[0] , vmax = z[1])
plt.colorbar()

#plt.savefig('centroids.pdf')
#%%   Plotting brightest pixel in each stamp vs wavlength

# Finding brightest pixel in each scut stamp
mpix = np.zeros(len(scut))
pixwl = np.zeros(len(scut))
indx = np.zeros(len(scut))
indy = np.zeros(len(scut))

for i in range(len(scut)):
    # Finding max value in each stamp
    mpix[i] = np.amax(scut[i])
    
    # Finding it's corresponding wavelength
    ind = np.unravel_index(np.argmax(scut[i], axis=None), scut[i].shape)
    ind = np.array((ind[0],ind[1]), dtype = int)
    pixwl[i] = wcut[i,ind[0],ind[1]]

# Subtracting background from it. I will use the bigger stamps (window) to do this.
# There is much more space around the star to get a good estimate of the background.
piece = np.zeros(len(scut))

for i in range(len(window)):
    piece[i] = np.average(scut[i,6:8,4:8])
background = np.average(piece)

#Plotting

plt.figure()
plt.scatter(pixwl, mpix - background, marker = '.')
plt.ylabel("Radiance [ergs/s/cm$^{2}$/$\AA$/sr]", fontsize = 12)
plt.xlabel("$\lambda$[$\mu$m]", fontsize = 12)
plt.title('Brightest Pixel vs. Wavelength')

#%% Finding omega_iso-lambda

R = 30 #resolution we want

files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Past Classes/RIT Undergrad/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data = image[0].data
waveext = image[1].data
wavedata = waveext[0,0:200,:]

iso = []
wv = np.zeros(13)

wv[0] = 1.7
for i in range(len(wv)-1):

    length = wv[i]/R
    wv[i+1] = wv[i] + length


diff = np.zeros((len(wv),len(wavedata)))
index_min= np.zeros((len(wv)), dtype = int)
minimum = np.zeros(len(wv))
index_min = np.zeros(len(wv))

for k in range(len(wv)):
    for j in range(len(wavedata)):
    
        diff[k,j] = abs(wv[k] - wavedata[j,0])
        
    index_min[k] = np.argmin(diff[k,:]) #stores indices

index_min = index_min[::-1].astype(int)

for h in range(len(index_min)-1):
    iso.append(wavedata[index_min[h]:index_min[h+1],:])


omega_iso = []
for r in iso:
    area = len(r)*len(r[0])
    omega_iso.append(area)

omega_iso = np.array([omega_iso])




#%%             Calculating flux in each box   

regrid_psf = regrid(epsf.data, 8,8)
# plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')

# Beam is normalized psf model source
beam = regrid_psf/np.nansum(regrid_psf)
outflux = np.zeros(len(scut))
wavelength = np.zeros(len(wcut))
err = np.zeros(len(errcut))
scut1 = scut
#wlwindow = wlwindow*10**4

for i in range(len(scut)):
    # Original data (no conversions)
    # outflux[i] = np.nansum(beam*scut1[i])/np.nansum(beam)
    # wavelength[i] = np.nansum(beam*wcut[i])/np.nansum(beam)
    # err[i] = np.nansum(beam*errcut[i])/np.nansum(beam)
    
    # Converted data
    outflux[i] = np.nansum(beam*scut1[i]*10**6)/np.nansum(beam)
    wavelength[i] = np.nansum(beam*wcut[i])/np.nansum(beam)
    outflux[i] = outflux[i]*wavelength[i]*10**4 # multiplying by lambda in A to get nW/m^2/sr
    err[i] = np.nansum(beam*errcut[i]*10**6)/np.nansum(beam)
    err[i] = err[i]*wavelength[i]*10**4

# diff = np.zeros((len(wv),len(outflux)))
# index_min= np.zeros((len(wv)), dtype = int)
# minimum = np.zeros(len(wv))
# index_min = np.zeros(len(wv))

# for k in range(len(wv)):
#     for j in range(len(wavelength)):
    
#         diff[k,j] = abs(wv[k] - wavelength[j])
        
#     index_min[k] = np.argmin(diff[k,:]) #stores indices

# index_min = index_min[::-1].astype(int)

# for h in range(len(index_min)-1):
#     err[index_min[h]:index_min[h+1]] = err[index_min[h]:index_min[h+1]]*errfactor[0,h]
 
    
plt.figure()
plt.errorbar(wavelength,outflux, yerr= err,capsize = 3,marker = '.',linestyle = 'None')
plt.ylabel("$\lambda I_{\lambda}$ [nW/m^2/sr]", fontsize = 12)
plt.xlabel("$\lambda$[$\mu$m]", fontsize = 12)
plt.ticklabel_format(axis="y", style="sci",scilimits=(-2,2))
plt.title('Vega Spectrum with LEISA', fontsize = 12)
#plt.savefig('LEISAerr.pdf')





#%%
import scipy.optimize as opt
def twoD_Gaussian( data_tuple, x0, y0, a, cx, cy, d):
    
    (x,y) = data_tuple
    x0 = float(x0) # x0 and y0 are b
    y0 = float(y0)
    
    g = d + a*np.exp( - (((x-x0)**2 / (2*cx**2))+((y-y0)**2/(2*cy**2))))
    return g.ravel()

xmin, xmax, nx = 0, 36, 37
ymin, ymax, ny = 0, 36, 37
x, y = np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)
p0 = [(20,18,0.08,0.5, 0.5, 0)]

data = epsf.data
popt, pcov = opt.curve_fit(twoD_Gaussian, (X,Y), data.ravel(), p0 = p0)
perr= np.sqrt(np.diag(pcov))
data_fitted = twoD_Gaussian((X,Y), *popt)
    
sigma_x = popt[3]
ds_x = perr[3]
    
sigma_y = popt[4]
ds_y = perr[4]
    
sigma_t = (sigma_x**2 + sigma_y**2)**0.5
ds_t = (((sigma_x*ds_x)**2+(sigma_y*ds_y)**2)**0.5)/np.sqrt(sigma_x**2 + sigma_y**2)
PSF_FWHM = 2*sigma_t*(2*np.log(2))**0.5
dPSF_FWHM = 2*ds_t*(2*np.log(2))**0.5



interval = ZScaleInterval() 
z = interval.get_limits(data)
fig, ax = plt.subplots(1, 1)

#plt.imshow(data,cmap='plasma',norm = None, vmin = z[0] , vmax = z[1])
plt.imshow(data_fitted.reshape(37,37), cmap='plasma', origin='upper',vmin = np.nanmin(data_fitted.reshape(37,37)) , vmax = np.nanmax(data_fitted.reshape(37,37)))
ax.contour(X, Y, data_fitted.reshape(37,37), 10, colors='black')
plt.colorbar()
#plt.savefig('Guassian.pdf')
plt.show()

fig = plt.figure(figsize=(20,10))
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y,data_fitted.reshape(37,37) ,cmap='viridis', edgecolor='none')
ax.contour(X, Y, data_fitted.reshape(37,37), 10, colors='black')
ax.view_init(25, -45)
plt.xlabel ('x');plt.ylabel ('y')
#plt.savefig('3dguassian.pdf')
plt.show()



#%%     Fitting Gaussians to stamps
import scipy.optimize as opt


def twoD_Gaussian( data_tuple, x0, y0, a, cx, cy, d):
    
    (x,y) = data_tuple
    x0 = float(x0) # x0 and y0 are b
    y0 = float(y0)
    
    g = d + a*np.exp( - (((x-x0)**2 / (2*cx**2))+((y-y0)**2/(2*cy**2))))
    return g.ravel()

xmin, xmax, nx = 0, 7, 8
ymin, ymax, ny = 0, 7, 8
x, y = np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)
p0 = [(4,4,0.003,0.5, 0.5, 0)]

sigma_x = np.zeros(len(scut))
ds_x = sigma_x*0

sigma_y = np.zeros(len(scut))
ds_y = sigma_y*0

sigma_t = sigma_x*0
ds_t = sigma_t*0
FWHM = sigma_x*0
dFWHM = FWHM*0

for i in range(len(scut)):
    data = scut[i]
    popt, pcov = opt.curve_fit(twoD_Gaussian, (X,Y), data.ravel(), p0 = p0)
    perr= np.sqrt(np.diag(pcov))
    data_fitted = twoD_Gaussian((X,Y), *popt)
    
    sigma_x[i] = popt[3]
    ds_x[i] = perr[3]
    
    sigma_y[i] = popt[4]
    ds_y[i] = perr[4]
    
    sigma_t[i] = (sigma_x[i]**2 + sigma_y[i]**2)**0.5
    ds_t[i] = (((sigma_x[i]*ds_x[i])**2+(sigma_y[i]*ds_y[i])**2)**0.5)/np.sqrt(sigma_x[i]**2 + sigma_y[i]**2)
    FWHM[i] = 2*sigma_t[i]*(2*np.log(2))**0.5
    dFWHM[i] = 2*ds_t[i]*(2*np.log(2))**0.5
    
x = np.linspace(0,len(scut),len(scut))
plt.figure()
plt.scatter(x, FWHM*60.83, marker = '.')
plt.xlabel('Stamp')
plt.ylabel('FWHM ($\mu$rad)')
plt.savefig('FWHMstamps.pdf')
plt.show()

# solid angle of the beam
pscale = 60.83e-6 # plate scale of LEISA rad/pixel
sang = 1.13*(FWHM*pscale)**2



interval = ZScaleInterval() 
z = interval.get_limits(data)
fig, ax = plt.subplots(1, 1)

#plt.imshow(data,cmap='plasma',norm = None, vmin = z[0] , vmax = z[1])
plt.imshow(data_fitted.reshape(8,8), cmap='plasma', origin='upper',vmin = np.nanmin(data_fitted.reshape(8,8)) , vmax = np.nanmax(data_fitted.reshape(8,8)))
ax.contour(X, Y, data_fitted.reshape(8,8), 10, colors='black')
plt.colorbar()
#plt.savefig('Guassian.pdf')
plt.show()

# 3-D plot

fig = plt.figure(figsize=(20,10))
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y,data_fitted.reshape(8,8) ,cmap='viridis', edgecolor='none')
ax.contour(X, Y, data_fitted.reshape(8,8), 10, colors='black')
ax.view_init(25, -45)
plt.xlabel ('x');plt.ylabel ('y')
#plt.savefig('3dguassian.pdf')
plt.show()

# sigma_x = popt[3]
# sigma_y = popt[4]
# sigma_t = (sigma_x**2 + sigma_y**2)**0.5
# FWHM = 2*sigma_t*(2*np.log(2))**0.5

# # solid angle of the beam
# pscale = 60.83e-6 # plate scale of LEISA rad/pixel
# sang = 1.13*(FWHM*pscale)**2



#%% Taking Kurucz reference and converting to nW/m^2/sr

vegafits =sorted(glob('C:/Users/Dennis Houlihan/Documents/Past Classes/RIT Undergrad/Capstone/Capstone 2/*.fits'))

image = fits.open(vegafits[0])
data = image[1].data
header = image[0].header




i = 1937
j = 2217
wl = data['WAVELENGTH'] 
wl = wl[i:j]
flux = data['FLUX']
flux = flux[i:j]
omega = 1.13*(np.average(FWHM)*60.83e-6)**2 # solid angle of the beam

omegap = 1.13*(np.average(PSF_FWHM))**2
domega = 2.26*np.average(FWHM)*np.average(dPSF_FWHM)

flux1 = flux*wl*10**(6)/omega
wl = wl*10**-4
plt.figure()
plt.plot(wl, flux1, color = 'black',label = 'Kurucz')
plt.xlabel('$\lambda$ ($\mu$m)')
plt.ylabel('$\lambda I_{\lambda}$ [nW/m^2/sr]')
#plt.scatter(wavelength,outflux, color = 'r', marker = '.', label = 'LEISA')
plt.legend(loc = 'best')
plt.savefig('Kmodelog.pdf')
plt.show()



#%% Fitting Kurucz to LEISA


store = np.zeros((len(wl),1))            # Creating empty arrays needed later
storeindex = np.zeros((len(wl),2)) 
minimum = np.zeros((len(wavelength),2), dtype = int)


for i in range(len(wavelength)):
    for j in range((len(wl))):

                    # Calculating distance from each object (Distance formula)
        store[j,0] = abs(wavelength[i] - wl[j])
        storeindex[j,1] = j #remembering index
        
   # minimum[i,0] = min(store)   #storing minimum values from store
    index_min1 = np.argmin(store) #stores indices
    minimum[i,0] = i #corresponds to the TEST/Bright Data indices
    minimum[i,1] = storeindex[index_min1,1] #corresponds to USNO DATA indices

wlnew = wavelength*0
flux1new = wavelength*0

for k in range(len(minimum)):
    wlnew[k] = wl[minimum[k,1]]
    flux1new[k] = flux1[minimum[k,1]]



def func(x,m,b):
     return m*x + b

#fit the data
#The sigma = means you will be doing a fit weighted by the uncertainty in y
par, cov = opt.curve_fit(func, flux1new, outflux)
perr= np.sqrt(np.diag(cov))
plt.figure()
plt.errorbar(flux1new,outflux, yerr= err, marker = '.', capsize = 3, linestyle = 'None')
plt.plot(flux1new,func(flux1new,*par), linewidth = 2,label = 'fit', color='black')
plt.ylabel('LEISA Data [nW/m$^{2}$/sr]')
plt.xlabel('Kurucz Model [nW/m$^{2}$/sr]')
#plt.savefig('LvsKline.pdf')
plt.show()
print('slope:', par[0])
print('offset:', par[1])
print('slope error:', perr[0])
print('Offset error:', perr[1])



#%% Taking difference between scaled model and data

outfluxerr = ((0.08e5)**2+(err)**2)**0.5
plt.figure()
plt.plot(wl, flux1*0.2, color = 'black',label = 'Kurucz', linewidth = 2.5)
plt.xlabel('$\lambda$ ($\mu$m)')
plt.ylabel('$\lambda$I$_{\lambda}$ (nW/m$^{2}$/sr)')
plt.errorbar(wavelength,outflux + 2.29e5, yerr= outfluxerr,capsize = 3,marker = '.',linestyle = 'None', label = 'LEISA Data')
#plt.scatter(wavelength,outflux+ 2.29e5, color = 'r', marker = '.', label = 'LEISA')
plt.legend(loc = 'best')
plt.title('Data after Scaling')

#%% Making Residual Plot

flux2new = flux1new*0.2
flux2newerr = flux2new*0.01

outfluxnew = outflux + 2.29e5
outfluxnewerr = ((0.08e5)**2+(err)**2)**0.5

difference = outfluxnew - flux2new
differr = ((flux2newerr)**2+(outfluxnewerr)**2)**0.5

RMS = (np.sum(difference**2)/(len(difference)))**0.5
RMSerr = np.sum(differr)*((0.5*2*np.sum(difference))/len(difference))*(np.sum(difference**2)/(len(difference)))**-0.5
print(RMS, RMSerr)


plt.figure()
plt.scatter(wavelength,difference, marker = '.')
#plt.errorbar(wavelength,difference, yerr = differr, marker = '.',capsize = 2, linestyle = 'None')
plt.hlines(0, wavelength[-1]-0.05,wavelength[0]+0.05, color =  'black')
plt.ylabel('Residuals from the Kurucz Model [nW/m$^{2}$/sr]')
plt.xlabel('$\lambda$ ($\mu$m)')
plt.xlim((1.68,2.38))

#%%    Binning the Residuals and finding sensitivity in each bin
flux2new = flux1new*0.2
flux2newerr = flux2new*0.01

outfluxnew = outflux + 2.29e5
outfluxnewerr = ((0.08e5)**2+(err)**2)**0.5

RMS = np.zeros(len(wv-1))
difference = []

diff = np.zeros((len(wv),len(wavelength)))
minimum = np.zeros(len(wv))
index_min = np.zeros(len(wv))

for k in range(len(wv)):
    for j in range(len(wavelength)):
    
        diff[k,j] = abs(wv[k] - wavelength[j])
        
    index_min[k] = np.argmin(diff[k,:]) #stores indices

index_min = index_min[::-1].astype(int)
poo = np.zeros(len(index_min))
#for i in range(len(RMS)):
for h in range(len(index_min)-1):
    difference.append(flux2new[index_min[h]:index_min[h+1]] - outfluxnew[index_min[h]:index_min[h+1]])
    poo[h] = wavelength[index_min[h]]
    
    
Std = np.zeros(len(difference))
N = np.zeros(len(difference))
for i in range(len(Std)):
    Std[i] = np.std(difference[i])
    N[i] = len(difference[i])

Std = Std[::-1]
N = N[::-1]
binned = Std*0
for i in range(len(wv)-1):
    binned[i] = (wv[i+1] + wv[i])/2 
    
t_int = 1.5
sensfactor = (omegap/omega_iso)**0.5
alpha = Std*np.sqrt(t_int/N)
dalpha = alpha/np.sqrt(N)

plt.figure()
plt.scatter(binned, alpha)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel('$\lambda [\mu m]$', fontsize = 12)
plt.ylabel('$\delta \lambda I_{\lambda}  [nW/m^{2}/sr/s^{1/2}]$', fontsize = 12)

# difference = flux2new - (outfluxnew)
# RMS = (np.sum(difference**2)/(len(difference)))**0.5


# Applying sqrt(omega_beam/omega_iso-lambda) factor to alpha
sensfactor = (omegap/omega_iso)**0.5
alphapri = alpha*sensfactor
alphapri_err = ((dalpha*np.sqrt(omegap/omega_iso))**2 + (alpha*domega/(2*np.sqrt(omega_iso*omegap)))**2)**0.5

plt.figure()
plt.errorbar(binned,alphapri[0,:], yerr= alphapri_err[0,:], marker = '.', capsize = 3, linestyle = 'None')
#plt.scatter(binned, alphapri)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel('$\lambda [\mu m]$', fontsize = 12)
plt.ylabel(r'$\alpha  [nW/m^{2}/sr/s^{1/2}]$', fontsize = 12)

#%%
import matplotlib.ticker as ticker
#from pylab import *
i = 1.69
j = wv[10]+0.008
scale_y = 1e5
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))




fig, axs = plt.subplots(2,1, figsize = (8,6), sharex=True)

axs[0].plot(wl, flux1*0.2, color = 'black',label = 'Kurucz', linewidth = 2.5)
axs[0].set_xlim(i,j)
axs[0].errorbar(wavelength,outfluxnew, yerr= outfluxerr,capsize = 3,marker = '.',linestyle = 'None', label = 'LEISA')
axs[0].vlines(wv, ymin = min(outflux), ymax = max(outflux), linestyle = '--', color = 'gray', label = 'Wavelength Bins')
#axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[0].yaxis.set_major_formatter(ticks_y)
#axs[0].text(0.0, 2.01, '1e5', fontsize=10, transform = gca().transAxes)
axs[0].axes.get_xaxis().set_visible(False)
axs[0].set_ylabel('$\lambda$I$_{\lambda}$ [10$^{5}$ nW m$^{-2}$ sr$^{-1}$]', fontsize = 11.5)

axs[1].errorbar(binned, alphapri[0,:], yerr = alphapri_err[0,:], marker = '.', linestyle = 'None', capsize = 3, color = 'black')
axs[1].set_xlim(i,j)
axs[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[1].yaxis.set_major_formatter(ticks_y)
axs[1].set_xlabel('$\lambda$ [$\mu$m]', fontsize = 11.5)
axs[1].set_ylabel('$\delta \lambda I_{\lambda}$ [10$^{5}$ nW m$^{-2}$ sr$^{-1}$ s$^{1/2}$]', fontsize = 11.5)
axs[0].legend(loc='best',prop={'size': 8})
plt.subplots_adjust(hspace=.0)
plt.savefig('AASpaperplot.pdf')







