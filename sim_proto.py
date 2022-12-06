# deboosting and sample completeness
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
from astropy.io import fits
import os
#import emcee
from astropy.convolution import Gaussian2DKernel
from astropy.modeling.models import Gaussian2D
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
from MF import *

# DZ: this script is to either evaluate the flux loss or generate 10,000 simulated maps. 
# Issue in terminal: python3 sim_proto.py
# it requires the rms map of the data and assumes you use the crop_map.py to generate the fits file, or rms map is the 2nd layer of your fits file. The radius of the map is hard coded and need to change the code directly

# The parameters for the fake source injection
large_s = 20
s_min = 0.05
step = 0.001

# generate a mock noise map
def mock_noise(sigma):
    '''returns the noise value with input sigma'''
    mu = 0
    mock = np.random.normal(mu,sigma)
    return mock

def mock_noisemap(imgg):
    '''return the noise map for a certain rms map'''
    mnp = np.empty(np.shape(imgg))
    for y in range(imgg.shape[0]):
        for x in range(imgg.shape[1]):
            mnp[y][x] = mock_noise(imgg[y][x]) # remember the 0 dim is the y axis and 1 dim is the x axis
    return mnp

def pdf_in(radius):
    '''number counts from Geach et al 2017, gives
    expected number of counts at each flux step, 
    hard coded a overdensity of 2 for bright 
    sources in the inner region'''
    r_deg = radius/60
    '''approx area'''
    area = math.pi*r_deg**2

    dn_norm = '(7180/2.5)*(s/2.5)**(-1.5)*(math.exp(-s/2.5))'

    n = [eval(dn_norm)*(step)*area for s in np.arange(s_min, large_s+step, step) if s <6]
    n1 = [3*eval(dn_norm)*(step)*area for s in np.arange(s_min, large_s+step, step) if s >=6]
    n = np.concatenate((n,n1))
    s = [s for s in np.arange(s_min, large_s+step, step)]

    a = np.array(n)
    summ = a.sum()
    return n, s, summ

def pdf_out(radius,radius_in):
    '''number counts from Geach et al 2017, gives
    expected number of counts at each flux step,
    hard coded a overdensity of 1 for bright 
    sources in the inner region'''
    r_deg = radius/60
    r_in = radius_in/60
    '''approx area'''
    area = math.pi*(r_deg**2-r_in**2)

    dn_norm = '(7180/2.5)*(s/2.5)**(-1.5)*(math.exp(-s/2.5))'

    n = [eval(dn_norm)*(step)*area for s in np.arange(s_min, large_s+step, step) if s <6]
    n1 = [2*eval(dn_norm)*(step)*area for s in np.arange(s_min, large_s+step, step) if s >=6]
    n = np.concatenate((n,n1))
    s = [s for s in np.arange(s_min, large_s+step, step)]

    a = np.array(n)
    summ = a.sum()
    return n, s, summ


def add_source(imgg, source, position, FWHM=12.6, scale=2):
    '''injected fake sources on the mock noise map'''
    '''
    fake_im = np.zeros_like(imgg)
    psf_data, lb1, lb2 = read_psf(FWHM=12.6, scale=2)
    x, y = position
    scale = source / psf_data[-lb1][-lb2]
    for i in range(psf_data.shape[0]):
        for j in range(psf_data.shape[1]):
            if 0<y+lb1+i<imgg.shape[0] and 0<x+lb2+j<imgg.shape[1]:
                fake_im[y+i+lb1][x+j+lb2]=scale*psf_data[i][j] # this is not very accurate, since the position of injected source is not purely random. I should make the position continuous rather than discrete and convolve it with the PSF.
    return fake_im
    '''
    sigma = FWHM / 2.355 / scale
    x, y = position
    gauss = Gaussian2D(source,x,y,sigma,sigma)
    gauss.render(imgg)

    

def insert(random_num):
    '''return the simulated map and associated catalog, radius is hard coded'''
    radius_in = 4
    radius = 9  # depend on the size of our map
    
    # Separate the NC for inner 4 and outer 5
    num_in, flux_in, summ_in = pdf_in(radius_in) 
    prob_in = [float(n / summ_in) for n in num_in]
    sources_in = np.random.choice(flux_in, size=int(summ_in), p=prob_in)
    number_in = len(sources_in) 
    
    num_out, flux_out, summ_out = pdf_out(radius,radius_in) 
    prob_out = [float(n / summ_out) for n in num_out]
    sources_out = np.random.choice(flux_out, size=int(summ_out), p=prob_out)
    number_out = len(sources_out) 
    
    img,header = crop(img0,header=head_ori,wcs=wcs_ori,radius=(radius+1)*60) #leave more space for the proper MF
    image_new = mock_noisemap(img)
    wcs = WCS(header)
    #noise_2 = image_new
    
    pri_spu, snr_spu, rms_spu, wcs_new_spu = smooth(file='./mock_850/spu_map'+str(random_num)+'.fits', img=image_new, noise=img, head=header,scale=2,FWHM=12.6,B_FWHM=26,radius=radius*60,mode=True) # the spurious map

    row = image_new.shape[0]
    col = image_new.shape[1]
    '''
    inj_in = []
    inj_out = []
    for i in range(row): # y axis
        for j in range(col): # x axis 
            if (i-int(row/2))**2+(j-int(col/2))**2 < (radius_in*30)**2: # due to the shape of the data 
                inj_in.append((j,i))
            elif (i-int(row/2))**2+(j-int(col/2))**2 < (radius*30)**2: # due to the shape of the data 
                inj_out.append((j,i))
    
    positions_in = random.sample(inj_in,int(number_in)) # source positions of inner region
    positions_out = random.sample(inj_out,int(number_out)) # outskirt source positions 
    
    # add sources together
    positions = np.concatenate((positions_in,positions_out))
    sources = np.concatenate((sources_in,sources_out))
    
    for i in range(len(sources)):
        image_new = image_new + add_source(sources[i],image_new,positions[i])
    '''
    x_pos = np.random.uniform(0,col,50000)
    y_pos = np.random.uniform(0,row,50000)
    positions = []
    positions.append(x_pos)
    positions.append(y_pos)
    positions = np.array(positions).reshape(50000,2)
    
    positions_in = np.array([po for po in positions if (po[0]-(col-1)/2)**2+(po[1]-(row-1)/2)**2 < (radius_in*30)**2])[:number_in] # po[0] is x axis and po[1] is y axis
    
    positions_out = np.array([po for po in positions if (radius_in*30)**2<=(po[0]-(col-1)/2)**2+(po[1]-(row-1)/2)**2 < (radius*30)**2])[:number_out]
    
    # add sources together
    positions = np.concatenate((positions_in,positions_out))
    sources = np.concatenate((sources_in,sources_out))
    
    #inj = Gaussian2D(0)
    
    for i in range(len(sources)):
        add_source(image_new,sources[i],positions[i])
        
        
    pri_mock, snr_mock, rms_mock, wcs_new = smooth(file='./mock_850/mock_map'+str(random_num)+'.fits', img=image_new, noise=img, head=header,scale=2,FWHM=12.6,B_FWHM=26,radius=radius*60,mode=True) # the mock map
    
    '''to record the inject, recovered and spurious sources'''
    err = []
    ra = []
    dec = []
    #snr0 = []
    err1 = []
    ra1 = []
    dec1 = []
    snr01 = []
    sources1 = []

    for fl, pos in zip(sources, positions):
        i, j = pos
        #err.append(rms_mock[i,j])
        ra.append(wcs.pixel_to_world(i,j).ra.degree/15)
        dec.append(wcs.pixel_to_world(i,j).dec.degree)
        apos = wcs_new.wcs_world2pix(wcs.pixel_to_world(i,j).ra.degree,wcs.pixel_to_world(i,j).dec.degree,1)
        xx = int(apos[0])
        yy = int(apos[1])
        err.append(rms_mock[yy,xx])   
        #snr0.append(fl/1000/img[i,j])
        if rms_mock[yy,xx] != 0: 
            if fl/rms_mock[yy,xx]> 3.5: # record the injected sources above the threshold in theory
                sources1.append(fl)
                err1.append(rms_mock[yy,xx])
                ra1.append(wcs.pixel_to_world(i,j).ra.degree/15)
                dec1.append(wcs.pixel_to_world(i,j).dec.degree)   
                snr01.append(fl/rms_mock[yy,xx]) # the supposed SNR
    # could also record two parts separately, but since only difference is the number of bright sources, and distance also kind of couple with noise, record them together. 
    t = Table([ra, dec, sources, err], names=['ra', 'dec', 'flux', 'err'])
    t.write('./mock_850/mock_map'+str(random_num)+'.dat', format='ascii', overwrite=True)
    t1 = Table([ra1, dec1, sources1, err1, snr01], names=['ra', 'dec', 'flux', 'err', 'snr'])
    t1.write('./mock_850/mock_map'+str(random_num)+'_detect.dat', format='ascii', overwrite=True)
    recovered_sources = find_sources('mock_850/mock_map'+str(random_num)+'_MF.fits', 3.5,FWHM=12.6,B_FWHM=26, scale=2, write=False)
    spurious_sources = find_sources('mock_850/spu_map'+str(random_num)+'_MF.fits', 3.5,FWHM=12.6,B_FWHM=26, scale=2, write=False)
    if recovered_sources.shape[0] == 0:
        f = open('./mock_850/mock_map'+str(random_num)+'_rec.dat','a+')
        f.write('No source has been detected.')
        f.close()
    else:
        rec_ra, rec_dec, rec_x, rec_y, rec_fluxes, rec_err, rec_snr = zip(*recovered_sources)
        t = Table([rec_ra, rec_dec, rec_fluxes, rec_err, rec_snr], names=('ra', 'dec', 'flux', 'rec_err', 'rec_snr'))
        t.write('./mock_850/mock_map'+str(random_num)+'_rec.dat', format='ascii.basic', overwrite=True)
    if spurious_sources.shape[0] == 0:
        g = open('./mock_850/spu_map'+str(random_num)+'.dat','a+')
        g.write('No fake source has been detected.')
        g.close()
    else:
        spu_ra, spu_dec, spu_x, spu_y, spu_fluxes, spu_err, spu_snr = zip(*spurious_sources)
        tt = Table([spu_ra, spu_dec, spu_fluxes, spu_err, spu_snr], names=('ra', 'dec', 'spu_flux', 'spu_err', 'spu_snr'))
        tt.write('./mock_850/spu_map'+str(random_num)+'.dat', format='ascii.basic', overwrite=True)
    os.system('rm ./mock_850/mock_map'+str(random_num)+'_MF.fits')
    os.system('rm ./mock_850/spu_map'+str(random_num)+'_MF.fits')

def correct(rms,flux):
    '''This function is to evaluate the flux loss during the Matched filtering, 
    it generates a dat file with the flux loss factor for each flux bin'''
    rms,header = crop(rms,header=head_ori,wcs=wcs_ori,radius=9*60)
    co_factor = []
    for i in range(100):
        mock = mock_noisemap(rms)
        row = int(noise.shape[0]/2)
        col = int(noise.shape[1]/2)
        
        add_source(mock,flux,(col,row))
        rec = smooth(file='correct.fits', img=mock, noise=rms, head=header, scale=2, FWHM=12.6, B_FWHM=26,radius=9*60,mode=True)[0]

        os.system('rm correct_MF.fits')
        flux_rec = np.max(rec[row-10:row+10,col-10:col+10])
        co_factor.append(flux/flux_rec)
    q = open('correct_850.dat','a+')
    q.write(str(flux)+'\t'+str(np.mean(co_factor))+'\t'+str(np.std(co_factor))+'\n')
    q.close()

if __name__ == '__main__': 
    # load the noise map
    filename = input('Please enter the file name of your cropped map (2nd layer should be the variance map): ') or '4C23_850_cal_crop.fits' # to avoid the problem on edge
    ori = fits.open(filename)
    img0 = np.sqrt(ori[1].data)# rms of the raw data
    head_ori = ori[1].header
    wcs_ori = WCS(head_ori)
    fcf = input('Estimate the correction factor? Yes or no: ')
    if fcf == 'Yes':
        for i in range(40):
            correct(img0,i/2+3)
    else:
        pool=Pool(os.cpu_count())
        pool.map(insert,np.arange(0,10000,1))    
#insert(0)

    
