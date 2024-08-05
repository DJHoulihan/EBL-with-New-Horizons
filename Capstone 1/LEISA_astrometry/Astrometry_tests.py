# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:27:05 2020

@author: -
"""
from glob import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from LEISA_Astrometry import RA_DEC
np.warnings.filterwarnings('ignore')

#%%

# Call the fits files 
files =sorted(glob('/Users/dennishoulihan/Documents/Capstone/Plutocruise-data/*.fit'))
#files = sorted(glob(os.getcwd()+'/Plutoencounter-data/*.fit'))
#files = sorted(glob('/Users/dennishoulihan/Documents/Capstone/Plutoencounter-data/*.fit'))

astrometry = RA_DEC(files[2])

ra = astrometry[0][:]
dec = astrometry[1][:]

openf = fits.open(files[5]) #opening a fitsfile    **
header = openf[0].header
boresight = np.array([header['SPCBRRA'],header['SPCBRDEC']])
quaternion = openf[7].data
#print(quaternion)

frames = np.linspace(1,len(astrometry[0]),len(astrometry[0]))
d_ra = boresight[0] - ra[:,127,127]
d_dec = boresight[1] - dec[:,127,127]

plt.figure()
plt.scatter(frames,d_ra, marker = 'o', label = 'Difference in RA')
plt.scatter(frames,d_dec, marker = 'o', label = 'Difference in Dec')
plt.plot([], [], ' ', label="{}".format(os.path.basename(files[5])))  # **
plt.xlabel('Image Frame')
plt.ylabel('Difference in degrees')
plt.title('Boresight Position vs. Middle of Each Exposure')
plt.legend(loc = 'best')
#plt.savefig('lsb703_test(1).pdf')

#%%
'''
plt.figure()

plt.xlabel('Image Frame')
plt.ylabel('Difference in degrees')
plt.title('Comparison of Boresight Position and Middle of Each Exposure')
plt.legend(loc = 'best')

#%%

#
#                   Testing Vega Coordinates
#

vega = np.loadtxt('vega1.txt', unpack = False)


#
#                       Difference in RA and Dec
#
plt.figure()
plt.scatter(vega[:,0], vega[:,3], marker = '.', color = 'red', label = 'Difference in RA')
plt.scatter(vega[:,0], vega[:,4], marker = '.', color = 'black', label = 'Difference in Dec')
plt.ylabel('Difference in Degrees')
plt.xlabel('Image Frame')
plt.plot([], [], ' ', label="{}".format(os.path.basename(files[0])))  # **
#plt.text(40,0.45,'lsb_0397097519_0x53c_sci.fit',bbox=dict(facecolor='gray', alpha=0.2))
plt.title("Vega's Accepted vs. Calculated J2000 Position")
plt.legend(loc = 'best')
#plt.savefig('Vega_test.pdf')

#%%
#
#               Total Distance between calculated vs accepted vega position
#
#plt.figure()
#plt.scatter(vega[:,0], vega[:,5], marker = 'o', color = 'red')
#plt.plot([], [], ' ', label="{}".format(os.path.basename(files[0])))  # **
#plt.xlabel('Image Frame')
#plt.legend(loc = 'best')
#plt.ylabel('Angular Distance (degrees)')
#plt.title("Comparison of Vega's Accepted vs. Calculated J2000 Position")
#plt.text(40,0.6,'lsb_0397097519_0x53c_sci.fit',bbox=dict(facecolor='gray', alpha=0.2))

#%%

#
#               Boresight Test for Vega Image
#
test = sorted(glob(os.getcwd()+'/Vega-test/*'))

vega_boresight = RA_DEC(test[0])

frames = np.linspace(1,len(vega_boresight[0]),len(vega_boresight[0]))
openf = fits.open(test[0]) #opening a fitsfile    **
header = openf[0].header
boresight = np.array([header['SPCBRRA'],header['SPCBRDEC']])
quaternion = openf[7].data
ra = vega_boresight[0][:]
dec = vega_boresight[1][:]
d_ra = boresight[0] - ra[:,127,127]
d_dec = boresight[1] - dec[:,127,127]

plt.figure()
plt.scatter(frames,d_ra, marker = '.', label = 'Difference in RA')
plt.scatter(frames,d_dec, marker = '.', label = 'Difference in Dec')
plt.plot([], [], ' ', label="{}".format(os.path.basename(test[0])))  # **
plt.xlabel('Image Frame')
plt.ylabel('Difference in degrees')
plt.title('Boresight Position vs. Middle of Each Exposure')
plt.legend(loc = 'best')
plt.savefig('Vegabor_test.pdf')



'''