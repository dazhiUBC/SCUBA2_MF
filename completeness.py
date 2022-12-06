import numpy as np
import statistics
import math
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import sys
import os
import time
from multiprocessing import Pool
# DZ: The script is adapted from Thomas's code. need to type "python3 completeness.py" on ur terminal and follow the commands. It can be changed to another version that can do the completeness analysis for different large kernels (so long as you have the simulated maps)
# PS: The positional accuracy should be high, so the threshold should be well below the instrumental FWHM 

#Change if you want more precise re-extract
if __name__ == '__main__': 
    err_max = float(input('The max of error: '))
    err_min = float(input('The min of error: '))

    flux_max = float(input('The max of flux: '))
    flux_min = float(input('The min of flux: '))

    path = input('The name of the directory (e.g. mock_850): ') + '/'
    os.system('mkdir '+path[5:-1])
    n_trials = int(input('The number of mock maps: '))
    #inik = int(input('The initial kernel size: '))
    #bigk = int(input('The number of different big kernels: '))

    #May have to change this based on how errors/fluxes distributed
    #same as correct_sources.py
    flux_bins = np.linspace(flux_min, flux_max, 26)
    error_bins = np.linspace(err_min, err_max, 26)

    flux_step = (flux_max-flux_min)/25
    error_step = (err_max-err_min)/25

    '''max arcsecond separation between input and extracted
    source to qualify as a match'''
    thresh = int(input('The researching radius: '))

def get_inputs(inp_file, rec_file):
  '''return the input and recover catalog for a simulated map'''
  t = Table.read(inp_file, format='ascii')
  t = t[t['err']>0.1]
  inp_ra, inp_dec, inp_flux, err = t['ra'], t['dec'], t['flux'], t['err']

  try:
    t2 = Table.read(rec_file, format='ascii') #ra dec flux rec_err rec_snr
    t2 = t2[t2['rec_err']>0.1]
    rec_ra, rec_dec, rec_flux, rec_err = t2['ra'], t2['dec'], t2['flux'], t2['rec_err']
  except:
    rec_ra = []
    rec_dec = []
    rec_flux = []
    rec_err = []
  inp_ra = [ir*15 for ir in inp_ra]
  rec_ra = [rr*15 for rr in rec_ra]
  return inp_ra, inp_dec, inp_flux, err, rec_ra, rec_dec, rec_flux, rec_err


def match(inp_file, rec_file, completeness=False):
  '''match the input source with the recovered sources and return the matched catalog'''
  inp_ra, inp_dec, inp_flux, inp_noise, rec_ra, rec_dec, rec_flux, rec_err = get_inputs(inp_file, rec_file)
  try:
      inp_ra = list(inp_ra)
      inp_dec = list(inp_dec)
      inp_flux = list(inp_flux)
      inp_noise = list(inp_noise)

      data = zip(rec_ra, rec_dec, rec_flux, rec_err)
      data2 = sorted(data, key=lambda tup: tup[2], reverse=True)
      rec_ra, rec_dec, rec_flux, rec_err = zip(*data2)
      rec_coord = SkyCoord(np.array(rec_ra)*u.degree, np.array(rec_dec)*u.degree)

      '''loops through outputs'''
      catalog = []
      for r_c, r_f, r_n in zip(rec_coord, rec_flux, rec_err):
        '''searches inputs to find match to recovered'''
        inp_coord = SkyCoord(np.array(inp_ra)*u.degree, np.array(inp_dec)*u.degree)
        d2d = r_c.separation(inp_coord).arcsecond

        matches = []
        ind = 0
        for d in d2d:
          if d < thresh:
            i_f = inp_flux[ind]
            i_noise = inp_noise[ind]
            matches.append((i_f, i_noise, ind))
          ind+=1

        if len(matches) > 0:
          match = max(matches, key = lambda t: t[0])
          catalog.append((match[0], r_n, r_f))

          '''remove the matched value (avoid doubles)'''
          index = match[2]
          del inp_ra[index]
          del inp_dec[index]
          del inp_flux[index]
          del inp_noise[index]

      if completeness:
        for i_f, i_n in zip(inp_flux, inp_noise):
          catalog.append((i_f, i_n, -99))
  except: 
    True
    
  return catalog

def compile_catalogs(n_trials, completeness=False):
  '''combines all catalogs for completeness or deboosting purpose'''
  master_catalog = []
  for i in range(0, n_trials):
    try:
        inp_file = path+'mock_map'+str(i)+'.dat'
        rec_file =path+'mock_map'+str(i)+'_rec.dat'
        if completeness:
          catalog = match(inp_file, rec_file, completeness=True)
        else:
          catalog = match(inp_file, rec_file, completeness=False)
        master_catalog.extend(catalog)
    except Exception:
        pass
  return master_catalog

def completeness(n_trials):
  '''completeness analysis and generate a dat file for completeness value'''

  master_catalog = compile_catalogs(n_trials,  completeness=True)

  bins = []
  for fb in flux_bins:
    for eb in error_bins:
      rec_bin = []
      for inp_f, err, rec_f in master_catalog:
        if fb<= inp_f < fb+flux_step and eb <= err < eb+error_step:
          rec_bin.append(rec_f)
      bins.append((fb, eb,rec_bin))

  t_counts = []
  counts = []
  for fb, eb, rb in bins:
    t_counts.append(len(rb)) # rb is the number of sources in the master catalog within the bin
    count = 0
    for r in rb:
      if r != -99: #-99 means no matched, !=-99 is matched sources
        count+=1 # count means the number of matched sources
    counts.append((fb,eb, count))

  input_fluxes, errors, n_recovered = zip(*counts)

  fraction = []
  for inp_flux, err, n_rec, tcount in zip(input_fluxes, errors, n_recovered, t_counts):
    '''if less than 10 fake sources went into an input flux/error bin,
    assume either 0% or 100% of sources recovered. this is just to make the
    completeness surface look nice.'''
    if tcount < 10:
      if float(inp_flux / err) < 3.5: #Again, change 3.5 to desired SNR detection level
        fraction.append((inp_flux, err, 0.00))
      if float(inp_flux / err) > 3.5: #Again, change 3.5 to desired SNR detection level
        fraction.append((inp_flux, err, 1.00))
    '''if more than 10 fake sources went into input flux/error bin, divide the
    number of sources recovered by sources originally added in that bin'''
    if tcount > 10:
      fraction.append((inp_flux, err, float(n_rec/tcount)))
    
  inp_flux, noise, frac = zip(*fraction)
  t = Table([inp_flux, noise, frac], names=('input_flux', 'noise', 'fraction'))
  t.write(path[5:]+'completeness_values.dat', format='ascii.basic')
  return fraction



def deboosting(n_trials):
  '''deboosting analysis and generate a dat file for deboosting value'''

  master_catalog = compile_catalogs(n_trials, completeness=False)
  '''gets ratio of input flux (from add.py) to recovered flux (after map
  processing and source extracting) for given noise and observed flux
  value. this is the deboosting value'''
  ratios = [(i_f, r_f, r_n, float(i_f/r_f)) for i_f, r_n, r_f in master_catalog]

  final_bins = []
  for fb in flux_bins:
    for eb in error_bins:
      bins = []
      for i_f, r_f, r_n, r in ratios:
        if fb <= r_f < fb+flux_step and eb <= r_n < eb+error_step:
        #Change 0.25 to match grid step defined at beginning
          bins.append(r)
      if len(bins) > 10:
        '''if at least 10 fake sources went into input flux/error bin, the deboosting
        measurement is deemed usable, otherwise the bin is ignored'''
        mean = statistics.mean(bins)
        median = statistics.median(bins)
        u75 = np.percentile(bins ,75)
        l25 = np.percentile(bins ,25)
        std = np.std(bins)
        final_bins.append((fb, eb, mean, median,std,u75,l25, len(bins)))

  rec_flux, noise, mean, med,std, u75, l25, num = zip(*final_bins)
  t = Table([rec_flux, noise, mean, med,std, u75, l25, num], names=('observed_flux', 'noise', 'deboost_value_mean', 'deboost_value_median','std','upper75', 'lower25', 'number'))
  t.write(path[5:]+'deboosting_values.dat',overwrite = True, format='ascii.basic')
  return master_catalog

def run():
  # for multiprocessing purpose only
  completeness(n_trials)
  deboosting(n_trials)
    
if __name__ == '__main__':
    run()
    #pool = Pool(20)    
    #pool.map(run,range(bigk))
