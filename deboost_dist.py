import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
from MF import *
from multiprocessing import Pool
import os
from astropy.coordinates import SkyCoord
import pandas as pd
#Now it only works for 850, because no many sources in 450 and 3mm, but we can still try to generate the number for them

source_file = input('Please put the catalog file (e.g. sources.dat): ') or 'sources_4C23_850_cal_crop_MF.dat'
deboosting_file = input('Please put the deboosting value file (e.g. deboosting_value.dat): ') or '850/deboosting_values.dat'
completeness_file = input('Please put the completeness value file (e.g. completeness_value.dat): ') or '850/completeness_values.dat'
#mode = input('Upper limit, lower limit or simple deboosting? (up, low, "enter")') or 'no'

#err_max = float(input('The max of error: ') or 2.5)
#err_min = float(input('The min of error: ') or 0.5)

#flux_max = float(input('The max of flux: ') or 21)
#flux_min = float(input('The min of flux: ') or 1)

#flux_bins = np.linspace(flux_min, flux_max, 26)
#error_bins = np.linspace(err_min, err_max, 26)

#flux_step = (flux_max-flux_min)/25
#error_step = (err_max-err_min)/25

# define the position of the HzRGs
rg = SkyCoord(ra=316.811766,dec=23.529172, unit= u.deg)

os.system('mkdir sim')
path = 'sim/'

def deboost(flux,mu,sigma):
    '''generate a deboosted flux through probability'''
    mock = np.random.normal(mu,sigma)
    return flux * mock

def get_sources(source_file):
  t = Table.read(source_file, format='ascii')
  ids, ra, dec, flux, noise, snr = t['obj_id'], t['ra'], t['dec'], t['flux'], t['err'], t['snr']
  sources = zip(ids, ra, dec, flux, noise, snr)
  raw_n = len(ids)
  return sources, raw_n


def get_values(deboosting_file):
  t = Table.read(deboosting_file, format='ascii')
  obs_flux, noise, sigma, value = t['observed_flux'], t['noise'], t['std'], t['deboost_value_mean']
  return obs_flux, noise, value, sigma

def get_values2(completeness_file):
  t = Table.read(completeness_file, format='ascii')
  inp_flux, noise, value = t['input_flux'], t['noise'], t['fraction']
  return inp_flux, noise, value


def deboost_sources(source_file,deboosting_file,num):
  sources, raw_num = get_sources(source_file)
  obs_flux, n, val, sigma = get_values(deboosting_file)

  deboosted_sources = []
  for obj_id, ra, dec, flux, noise, snr in sources:
    ind1 = np.searchsorted(obs_flux, flux, side = 'left') # caution: this is i-1<v<i
    selected_flux = obs_flux[ind1-1]

    l_err, l_vals, l_sig = [], [], []
    for of, e, v, si in zip(obs_flux, n, val, sigma):
      if of == selected_flux:
        l_err.append(e)
        l_vals.append(v)
        l_sig.append(si)

    ind2 = np.searchsorted(l_err, noise, side = 'left')
    selected_noise = l_err[ind2-1]
    deboosting_value = l_vals[ind2-1]
    unc = l_sig[ind2-1]
    fl = deboost(flux,deboosting_value,unc)
    
    obj_id = str(obj_id).zfill(3)
    deboosted_sources.append((obj_id, ra, dec, flux, fl, noise))

  ids, ra, dec, flux, fl, noise  = zip(*deboosted_sources)
  t = Table([ids, ra, dec, flux, fl , noise ], names=('obj_id', 'ra', 'dec', 'origin_flux', 'deboost_flux', 'err' ))
  t.write(path+str(num)+'deboosted_'+source_file, format='ascii.basic')

  return str(str(num)+'deboosted_'+source_file)

def get_sources1(source_file):
  t = Table.read(source_file, format='ascii')
  ids, ra, dec, flux, dflux, noise = t['obj_id'], t['ra'], t['dec'], t['origin_flux'], t['deboost_flux'], t['err'] 
  sources = zip(ids, ra, dec, flux, dflux, noise)
  raw_n = len(ids)
  return sources, raw_n

def completeness_correct(flux,noise):
  inp_flux,noi , val = get_values2(completeness_file)
  #f_st = inp_flux[1]-inp_flux[0]
  completeness_value = 1
  for f,n, va in zip(inp_flux,noi,val):
    if f>=flux:
      if n>=noise:
        completeness_value = va
  n_val = 1/completeness_value
  return n_val

def correct_deboosted_sample(file):
  deb_sources, raw_n = get_sources1(path+file)
  ids, ra, dec, flux, dflux, noise = zip(*deb_sources)
  #total_sources = 0
  '''for deboosted flux and noise level, get number of sources
  we should have found if we didn't miss any(completeness correct)'''
  total_sources = np.zeros(15)
  inner_sources = np.zeros(15)
  for f, n, ra, dec in zip(dflux, noise, ra, dec):
    pos = SkyCoord(ra=ra*15,dec=dec, unit= u.deg)
    n_val = completeness_correct(f,n)
    flux_bi = []
    for i in range(15):
      flux_bi.append(i+1)
      if i+1<=f<i+2:
        total_sources[i] += n_val
        if pos.separation(rg).arcmin <=4:
          inner_sources[i] += n_val
    
  t = Table([flux_bi,total_sources,inner_sources], names=('flux', 'num', 'num_in'))
  t.write(path+'NC_'+file, format='ascii.basic')
  return 'NC_'+file

def fid_cor850(fname,fidname='fide_flux_850.npy'):
    # only correct through flux, if consider noise it might be biased (low flux & large noise never corrected)
    fid_flux = np.load(fidname)
    obs_flux = fid_flux[0]
    fidelity = fid_flux[1]
    source = pd.read_csv(path+fname,delimiter= ' ')
    co = np.zeros(len(source['flux']))
    co_in = np.zeros(len(source['flux']))
    if source['flux'][0] == obs_flux[0]:
        for i in range(len(source['flux'])):
            co[i] = source['num'][i] * fidelity[i]
            co_in[i] = source['num_in'][i] * fidelity[i]
        t = Table([source['flux'],co,co_in], names=('flux', 'cor_num','cor_num_in'))
        t.write(path+'Ture_'+fname, format='ascii.basic',overwrite=True)
    else:
        print('Please check the data')
        
def simulation(number): 
  '''deboost should be enoguh, completeness and fidelity correction should be done after the NC count'''  
  file = deboost_sources(source_file,deboosting_file,number)
  cfile = correct_deboosted_sample(file)
  fid_cor850(cfile)
'''
def nc(file):
  #havn't finished yet
  deb_sources, raw_n = get_sources1(path+file)
  ids, ra, dec, flux, dflux, noise = zip(*deb_sources)
  total_sources = 0
  #for deboosted flux and noise level, get number of sources
  #we should have found if we didn't miss any(completeness correct)
  total_sources = np.zeros(len(flux_bins))
  for f, n in zip(dflux, noise):
    n_val = completeness_correct(f, n)
    for i in range(len(flux_bins)):
      if flux_bins[i]<=f<flux_bins[i]+flux_step:
        total_sources[i] += n_val
    

  t = Table([flux_bins,total_sources], names=('flux', 'num'))
  t.write(path+'NC_'+file, format='ascii.basic')
  return 'NC_'+file
'''
def nc_reader(num):
    '''After obtain the individual NC, we use this fun to get the total NC '''
    file = 'sim/Ture_NC_'+str(num)+'deboosted_sources_4C23_850_cal_crop_MF.dat' 
    nc = pd.read_csv(file,delimiter=' ')
    return nc['cor_num']

if __name__=='__main__':
    pool=Pool(os.cpu_count())
    pool.map(simulation,range(10000))    
