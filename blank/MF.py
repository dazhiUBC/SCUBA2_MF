import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
from astropy.io import fits
import os
from astropy.convolution import Gaussian2DKernel
from multiprocessing import Pool
import math
import random
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import math
import datetime as DT
from astropy.table import Table
import sys
import os
from astropy.nddata import Cutout2D
from astropy.convolution import convolve
from astropy.wcs.utils import pixel_to_skycoord
from astropy import units as u

# DZ: To use this script for the matched filtering process, 
# issue in terminal: python3 MF.py
# It will ask you for your data, noise map and other relevant information

def crop(img,header,wcs,radius):
    '''crop the image to remove some junk, return the array and header'''
    side = radius*2.
    size = u.Quantity((side+10, side+10), u.arcsec) # make the image slightly larger
    position = pixel_to_skycoord(header['CRPIX1'], header['CRPIX2'], wcs)
    # Make the cutout, including the WCS
    cutout = Cutout2D(img, position=position, size=size, wcs=wcs)
    return cutout.data,cutout.wcs.to_header()#*(r<=radius/2)

def r_theta(im, xc, yc):
    '''make a radius mask, and returns the radius rr 
    and the angle phi for point (xc,yc)'''
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    return(rr, phi)

def bigpsf(FWHM = 30,scale=2):
    '''The size of the Big Gaussian kernel for background subtraction, 
    need to be sufficiently large, return psf, and the x and y size'''
    sigma = FWHM / 2.355 / scale # convert the FWHM to sigma
    Kernel = Gaussian2DKernel(sigma,x_size=51,y_size=51) # can either use to smooth the image or treat as psf
    psf_data = np.array(Kernel)
    lb1 = -int(Kernel.shape[0]/2) # The y axis  
    lb2 = -int(Kernel.shape[1]/2) # The x axis 
    return psf_data, lb1, lb2

def read_psf(FWHM = 9,scale=2):
    '''Return the peak normalized instrumentation PSF and x, y size. 
    If FWHM/scale is too large, kernel size might need to be changed'''
    sigma = FWHM / 2.355 / scale # convert the FWHM to sigma
    Kernel = Gaussian2DKernel(sigma,x_size=51,y_size=51) # can either use to smooth the image or treat as psf
    psf_data = np.array(Kernel)
    lb1 = -int(Kernel.shape[0]/2) # The y axis  
    lb2 = -int(Kernel.shape[1]/2) # The x axis 
    return psf_data/np.max(psf_data), lb1, lb2

def smooth(file,img,noise,head,scale,FWHM,B_FWHM,radius,mode=False):
    '''The whole MF process, removal of the large scale noise and do the 
    regular MF, returns the fits data. The first layer is science map,
    2nd is SNR and the last one is the noise map.'''
    wcs = WCS(head)
    P = read_psf(FWHM,scale)[0]
    P = P/np.max(P) # The instrumental PSF should be peak normalized
    B = bigpsf(B_FWHM,scale)[0] # The smoothing PSF should be volume normalized 
    
    '''smooth both map and psf with big gaussian
    Maintain the original rms since the smoothed rms is negligible'''
    Bimg = convolve(img,B,nan_treatment='fill') # it should be volume normalized,normalize_kernel=False)
    #Bvar = convolve(var,B**2,normalize_kernel=False)
    BP = convolve(P,B,nan_treatment='fill')#,normalize_kernel=False)
    
    '''Subtract the smoothed map'''
    img = img-Bimg
    P = P-BP 
        
    '''The normal MF process'''
    W = np.empty(np.shape(img))
    Wi = np.empty(np.shape(img))
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if noise[i][j]!=0:
                W[i][j]=1/(noise[i][j])**2
                Wi[i][j]=img[i][j]*W[i][j]
    F1 = convolve(Wi,P,normalize_kernel=False,nan_treatment='fill')
    F2 = convolve(W,P**2,normalize_kernel=False,nan_treatment='fill')
    F = F1/F2
    F3 = 1/np.sqrt(F2)
    FF = F1*F3
    
    '''Store the MF map with proper cropping'''
    r,t = r_theta(FF,int(head['CRPIX1']),int(head['CRPIX2']))
    pri,head_pri = crop(F* (r <= radius/scale),header=head,wcs=wcs,radius=radius)
    primaryhdu = fits.PrimaryHDU(pri, header = head_pri) 
    snr, head_snr = crop(FF* (r <= radius/scale),header=head,wcs=wcs,radius=radius)
    snrhdu = fits.ImageHDU(snr, header = head_snr)
    rms, head_rms = crop(F3* (r <= radius/scale),header=head,wcs=wcs,radius=radius)
    rmshdu = fits.ImageHDU(rms, header = head_rms)
    greyHDU=fits.HDUList([primaryhdu,snrhdu,rmshdu]) 
    greyHDU.writeto(file[:-5]+'_MF.fits')
    if mode:
        return pri, snr, rms, WCS(head_rms)

def find_sources(fileid, lim, FWHM=12.6,B_FWHM=26, scale=2, write=True):
    '''The source extraction algorithm is adapted originated from Thomas's script, returns the raw 
    catalog and save the dat file'''
    snr_map = fits.open(fileid)[1].data
    flux_map = fits.open(fileid)[0].data
    psf_data, lb1, lb2 = read_psf(FWHM,scale)
    bk_data = bigpsf(B_FWHM,scale)[0]
    BP = convolve(psf_data,bk_data,nan_treatment='fill')
    P = psf_data-BP 
    psf_data = P/np.max(P)
    #psf_data, lb1, lb2 = read_psf(FWHM,scale) PSF should be the subtracted PSF
    error_map = fits.open(fileid)[2].data 
    wcs = WCS(fits.open(fileid)[0].header)
  

    sources = []
    positions = []

    while np.nanmax(snr_map) > float(lim):
        brightest_row, brightest_col = np.where(np.nanmax(snr_map) == snr_map)
        brightest_row, brightest_col = int(brightest_row[0]), int(brightest_col[0])

        snr_val = snr_map[brightest_row, brightest_col]
    #print(snr_val)
        flux_val = flux_map[brightest_row, brightest_col]
        e_val = error_map[brightest_row, brightest_col]
    #print(e_val)

        snr_scale = snr_val / psf_data[lb1-1,lb2-1]
        flux_scale = flux_val / psf_data[lb1-1,lb2-1]

    #Below may have to be changed depending on map size. change if error occurs.
    #The error will be: ValueError: operands could not be broadcast together with shapes
    #(XXX, XXX) (XXX, XXX). Add or subtract an integer in blanks so they match maps
    #this occurs bc of imperfect cropping at edges.
        snr_blank = np.zeros_like(snr_map)
        flux_blank = np.zeros_like(flux_map)


        ind = []
        for i in range(-8,9):
            for j in range(-8,9):
                ind.append((i,j))

        for i,j in ind:
            try:
                snr_blank[brightest_row+i, brightest_col+j] = snr_scale*psf_data[lb1-1+i,lb2-1+j]
                flux_blank[brightest_row+i, brightest_col+j] = flux_scale*psf_data[lb1-1+i,lb2-1+j]
            except IndexError:
                continue

        snr_map = np.subtract(snr_map, snr_blank)
        flux_map = np.subtract(flux_map, flux_blank)

    # Get ra, dec from pixel positions
    # TRG: added +1 to _col and _row.
        wx, wy = wcs.wcs_pix2world(brightest_col+1, brightest_row+1, 1)

        sources.append([wx/15., wy, brightest_col+1, brightest_row+1, flux_val, e_val, snr_val])
    # TRG: added +1 to _col and _row.
        positions.append([brightest_col+1, brightest_row+1])

    sources = np.array(sources).reshape(np.array(sources).shape[:2])
    positions = np.array(positions).reshape(np.array(positions).shape[:2])
    p = Table(positions)

    if write:
        ids = [str(x).zfill(3) for x in range(1, len(sources[:])+1)]
        ra, dec, x, y, flux, flux_err, snr = zip(*sources)
        t = Table([ids, ra, dec, x, y, flux, flux_err, snr], names=['obj_id', 'ra', 'dec', 'x', 'y', 'flux', 'err', 'snr'])
        t.write('sources_'+fileid[:-5]+'.dat', format='ascii', overwrite=True)

    return sources

if __name__ == '__main__':  
    file = input('Please enter the file name of your original map: ')
    hdu = fits.open(file)
    img = (hdu[0].data)
    head = hdu[0].header
    nfile = input('Please enter the noise file name here, if the variance, weight, or rms are already in the second layer of the map, please specified:(i.e. variance, weight, rms) ')
    if nfile == 'variance':
        noise = np.sqrt(hdu[1].data)
    elif nfile == 'weight':
        noise = 1/np.sqrt(hdu[1].data)
    elif nfile == 'rms':
        noise = (hdu[1].data)
    else:
        noise = fits.open(nfile)[0].data

    scale = int(input('The pixel scale of the map is: ') )
    FWHM = float(input('The instrumental FWHM is: '))
    B_FWHM = int(input('The big Gaussian FWHM is: '))
    radius = int(input('Please enter the radius you want to maintain (in arcsec): '))   

    smooth(file=file,img=img,noise=noise,head=head,scale=scale,FWHM=FWHM,B_FWHM=B_FWHM,radius=radius)

    find_sources(file[:-5]+'_MF.fits',3.5,FWHM=FWHM,B_FWHM=B_FWHM, scale=scale)