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
from astropy.io import fits

# this is Dazh's new completeness code. The treatment for the spurious source estimation and effective area is corrected. 

def get_inputs(inp_file, rec_file):
    '''
    Return the input and recover catalog for a simulated map.
    '''
    t = Table.read(inp_file, format='ascii')
    t = t[(t['err'] > 0.1) & (t['flux'] > t['err'])]
    inp_ra = [ir * 15 for ir in t['ra']]
    inp_dec, inp_flux, err = t['dec'], t['flux'], t['err']

    rec_ra, rec_dec, rec_flux, rec_err = [], [], [], []
    try:
        t2 = Table.read(rec_file, format='ascii')
        t2 = t2[t2['rec_err'] > 0.1]
        rec_ra = [rr * 15 for rr in t2['ra']]
        rec_dec, rec_flux, rec_err = t2['dec'], t2['flux'], t2['rec_err']
    except:
        print('error in get_inputs')

    return inp_ra, inp_dec, inp_flux, err, rec_ra, rec_dec, rec_flux, rec_err

def match(inp_file, rec_file ):
    '''Match the input source with the recovered sources and return the matched, lost, spurious
    catalogs, with shape of (n,3), (n,2), and (n,2), respectively'''
    inp_ra, inp_dec, inp_flux, inp_noise, rec_ra, rec_dec, rec_flux, rec_err = get_inputs(inp_file, rec_file)
    matches = []
    lost = []
    spurious = []
    try:
        # Convert input lists to array
        inp_ra = np.array(inp_ra)
        inp_dec = np.array(inp_dec)
        inp_flux = np.array(inp_flux)
        inp_noise = np.array(inp_noise)
        inp_coord = SkyCoord( inp_ra * u.degree, inp_dec * u.degree)
        
        # Sort recovered sources by flux in descending order
        data = zip(rec_ra, rec_dec, rec_flux, rec_err)
        data2 = sorted(data, key=lambda tup: tup[2], reverse=True)
        rec_ra, rec_dec, rec_flux, rec_err = zip(*data2)
        rec_ra = np.array(rec_ra)
        rec_dec = np.array(rec_dec)
        rec_flux = np.array(rec_flux)
        rec_err = np.array(rec_err)
        rec_coord = SkyCoord( rec_ra * u.degree, rec_dec * u.degree)

        # Loop through recovered sources
        for i in range(len(rec_coord)):
            # Search inputs to find a match to the recovered source
            d2d = rec_coord[i].separation(inp_coord).arcsecond
            n_match = len(inp_coord[d2d<=thresh])
            
            if n_match > 0:
                if n_match == 1:
                    index = int(np.where(d2d<=thresh)[0])
                
                else: 
                    n_index = np.where(d2d<=thresh)[0]
                    can_flux = inp_flux[n_index]
                    index = n_index[np.argmax(can_flux)]
                    
                matches.append((inp_flux[index], rec_flux[i], rec_err[i]))
                    
                # Remove the matched value to avoid duplicates
                inp_ra = np.delete(inp_ra,index)
                inp_dec = np.delete(inp_dec,index)
                inp_flux = np.delete(inp_flux,index)
                inp_noise = np.delete(inp_noise,index)
                inp_coord = np.delete(inp_coord,index)
            elif n_match == 0:
                spurious.append((rec_flux[i],rec_err[i]))
        
        # if no counterpart is found in recovered sources, it means the source is lost
        for i in range(len(inp_flux)):
            lost.append((inp_flux[i],inp_noise[i])) 
    except:
        print('error in match')

    return matches, lost, spurious

def compile_catalogs(n):
    '''return the big match, lost, and spurious catalog and also the number of spurious sources'''
    master_m = []
    master_l = []
    master_s = []
    n_spu = []
    for i in range(n):
        try:
            inp_file = path+'mock_map'+str(i)+'.dat'
            rec_file =path+'mock_map'+str(i)+'_rec.dat'
            matches, lost, spurious = match(inp_file, rec_file)
            master_m.extend(matches)
            master_l.extend(lost)
            master_s.extend(spurious)
            n_spu.append(len(spurious))
        except:
            pass
    return np.array(master_m), np.array(master_l), np.array(master_s), n_spu

def analysis(n):
    '''perform completeness, deboosting, and fidelity analysis'''
    master_m, master_l, master_s, n_spu = compile_catalogs(n)
    completeness(master_m,master_l)
    deboosting(master_m)
    spurious(master_m,master_s,n_spu)
    return 0
    
def completeness(master_m,master_l):
    '''completeness analysis and generate a dat file for completeness value'''
    match_i = master_m[:,0]
    match_n = master_m[:,2]
    lost_i = master_l[:,0]
    lost_n = master_l[:,1]
    bins = []
    for fb in flux_bins:
        for eb in error_bins:
            n_rec = len(match_i[(match_i>=fb)*(match_i<fb+flux_step)*(match_n>=eb)*(match_n<eb+error_step)])
            n_los = len(lost_i[(lost_i>=fb)*(lost_i<fb+flux_step)*(lost_n>=eb)*(lost_n<eb+error_step)])
            n_sum = n_rec+n_los
            if n_sum < 10: 
                if fb/eb >=3.5:
                    frac = 1.
                else: 
                    frac = 0.
            else: 
                frac = n_rec/n_sum
            bins.append((fb,eb,frac,n_sum))
    inp_flux, noise, frac,num = zip(*bins)
            
    t = Table([inp_flux, noise, frac,num], names=('input_flux', 'noise', 'fraction', 'number'))
    t.write(path[8:]+'completeness_values.dat', format='ascii.basic')
    
    f_bins = []
    for fb in flux_bins:
        n_rec = len(match_i[(match_i>=fb)*(match_i<fb+flux_step)])
        n_los = len(lost_i[(lost_i>=fb)*(lost_i<fb+flux_step)])
        n_sum = n_rec+n_los
        if n_sum < 10: 
            if fb >=5:
                frac = 1.
            else: 
                frac = 0.
        else: 
            frac = n_rec/n_sum
        f_bins.append((fb,frac,n_sum))
    inp_flux, frac, num = zip(*f_bins)
            
    t = Table([inp_flux, frac, num], names=('input_flux', 'fraction', 'number'))
    t.write(path[8:]+'completeness_flux.dat', format='ascii.basic')
    
    e_bins = []
    for eb in error_bins:
        n_rec = len(match_i[(match_n>=eb)*(match_n<eb+error_step)])
        n_los = len(lost_i[(lost_i>=3.5*lost_n)*(lost_n>=eb)*(lost_n<eb+error_step)])
        n_sum = n_rec+n_los
        if n_sum < 10: 
            if eb <= 1:
                frac = 1.
            else: 
                frac = 0.
        else: 
            frac = n_rec/n_sum
        e_bins.append((eb,frac,n_sum))
    error, frac,num = zip(*e_bins)
            
    t = Table([error, frac,num], names=('noise', 'fraction', 'number'))
    t.write(path[8:]+'completeness_noise.dat', format='ascii.basic')
    
    snr_bins = []
    for i in range(40):
        n_rec = len(match_i[(match_i/match_n>=2.75+0.5*i)*(match_i/match_n<3.25+0.5*i)])
        n_los = len(lost_i[(lost_i/lost_n>=2.75+0.5*i)*(lost_i/lost_n<3.25+0.5*i)])
        n_sum = n_rec+n_los
        if n_sum < 10: 
            if i >= 5:
                frac = 1.
            else: 
                frac = 0.
        else: 
            frac = n_rec/n_sum
        snr_bins.append((0.5*i+3,frac,n_sum))
    inp_snr, frac, num = zip(*snr_bins)
            
    t = Table([inp_snr, frac, num], names=('input_snr', 'fraction', 'number'))
    t.write(path[8:]+'completeness_snr.dat', format='ascii.basic')
    return 0

def deboosting(master_m):
    '''deboosting analysis and generate a dat file for deboosting value'''
    match_i = master_m[:,0]
    match_r = master_m[:,1]
    match_n = master_m[:,2]
    ratios = match_i/match_r
    bins = []
    for fb in flux_bins:
        for eb in error_bins:
            deb = ratios[(match_r>=fb)*(match_r<fb+flux_step)*(match_n>=eb)*(match_n<eb+error_step)]
            n_sum = len(deb)
            if n_sum >= 10: 
                mean = np.mean(deb)
                med = np.mean(deb)
                u84 = np.percentile(deb,84)
                l16 = np.percentile(deb,16)
                std = np.std(deb)
                bins.append((fb, eb, mean, med, std, u84, l16, n_sum))
                
    rec_flux, noise, mean, med,std, u84, l16, num = zip(*bins)
    t = Table([rec_flux, noise, mean, med,std, u84, l16, num], names=('observed_flux', 'noise', 'deboost_value_mean', 'deboost_value_median','std','upper84', 'lower16', 'number'))
    t.write(path[8:]+'deboosting_values.dat',overwrite = True, format='ascii.basic')
    
    f_bins = []
    for fb in flux_bins:
        deb = ratios[(match_r>=fb)*(match_r<fb+flux_step)]
        n_sum = len(deb)
        if n_sum >= 10: 
            mean = np.mean(deb)
            med = np.mean(deb)
            u84 = np.percentile(deb,84)
            l16 = np.percentile(deb,16)
            std = np.std(deb)
            f_bins.append((fb, mean, med, std, u84, l16, n_sum))
                
    rec_flux, mean, med,std, u84, l16, num = zip(*f_bins)
    t = Table([rec_flux, mean, med,std, u84, l16, num], names=('observed_flux', 'deboost_value_mean', 'deboost_value_median','std','upper84', 'lower16', 'number'))
    t.write(path[8:]+'deboosting_flux.dat',overwrite = True, format='ascii.basic')
    
    e_bins = []
    for eb in error_bins:
        deb = ratios[(match_n>=eb)*(match_n<eb+error_step)]
        n_sum = len(deb)
        if n_sum >= 10: 
            mean = np.mean(deb)
            med = np.mean(deb)
            u84 = np.percentile(deb,84)
            l16 = np.percentile(deb,16)
            std = np.std(deb)
            e_bins.append((eb, mean, med, std, u84, l16, n_sum))
                
    noise, mean, med,std, u84, l16, num = zip(*e_bins)
    t = Table([noise, mean, med,std, u84, l16, num], names=('noise', 'deboost_value_mean', 'deboost_value_median','std','upper84', 'lower16', 'number'))
    t.write(path[8:]+'deboosting_noise.dat',overwrite = True, format='ascii.basic')
    
    snr_bins = []
    for i in range(40):
        deb = ratios[(match_r/match_n>=2.75+0.5*i)*(match_r/match_n<3.25+0.5*i)]
        n_sum = len(deb)
        if n_sum >= 10: 
            mean = np.mean(deb)
            med = np.mean(deb)
            u84 = np.percentile(deb,84)
            l16 = np.percentile(deb,16)
            std = np.std(deb)
            snr_bins.append((0.5*i+3, mean, med, std, u84, l16, n_sum))
                
    snr, mean, med,std, u84, l16, num = zip(*snr_bins)
    t = Table([snr, mean, med,std, u84, l16, num], names=('observed_snr', 'deboost_value_mean', 'deboost_value_median','std','upper84', 'lower16', 'number'))
    t.write(path[8:]+'deboosting_snr.dat',overwrite = True, format='ascii.basic')
    return 0

def spurious(master_m,master_s,n_spu):
    '''spurious analysis and generate a dat file for contamination value,
    Note: in the old code, I made a mistake. The fidelity correction should be done
    before the deboosting.'''
    match_r = master_m[:,1]
    match_n = master_m[:,2]
    fake_r = master_s[:,0]
    fake_n = master_s[:,1]

    bins = []
    for fb in flux_bins:
        for eb in error_bins:
            n_rec = len(match_r[(match_r>=fb)*(match_r<fb+flux_step)*(match_n>=eb)*(match_n<eb+error_step)])
            n_fake = len(fake_r[(fake_r>=fb)*(fake_r<fb+flux_step)*(fake_n>=eb)*(fake_n<eb+error_step)])
            n_sum = n_rec+n_fake
            if n_sum < 10: 
                if fb/eb >=5:
                    frac = 1.
                else: 
                    frac = 0.
            else: 
                frac = n_rec/n_sum
            bins.append((fb,eb,frac,n_sum))
    rec_flux, noise, frac, num = zip(*bins)

    t = Table([rec_flux, noise, frac, num], names=('rec_flux', 'noise', 'fidelity', 'number'))
    t.write(path[8:]+'fidelity_values.dat', format='ascii.basic')

    f_bins = []
    for fb in flux_bins:
        n_rec = len(match_r[(match_r>=fb)*(match_r<fb+flux_step)])
        n_fake = len(fake_r[(fake_r>=fb)*(fake_r<fb+flux_step)])
        n_sum = n_rec+n_fake
        if n_sum < 10: 
            if fb >=5:
                frac = 1.
            else: 
                frac = 0.
        else: 
            frac = n_rec/n_sum
        f_bins.append((fb,frac,n_sum))
    rec_flux,frac, num = zip(*f_bins)

    t = Table([rec_flux, frac, num], names=('rec_flux', 'fidelity', 'number'))
    t.write(path[8:]+'fidelity_flux.dat', format='ascii.basic')

    e_bins = []
    for eb in error_bins:
        n_rec = len(match_r[(match_n>=eb)*(match_n<eb+error_step)])
        n_fake = len(fake_r[(fake_n>=eb)*(fake_n<eb+error_step)])
        n_sum = n_rec+n_fake
        if n_sum < 10: 
            if eb >=1:
                frac = 0.
            else: 
                frac = 1.
        else: 
            frac = n_rec/n_sum
        e_bins.append((eb,frac,n_sum))
    noise,frac, num = zip(*e_bins)

    t = Table([noise, frac, num], names=('noise', 'fidelity', 'number'))
    t.write(path[8:]+'fidelity_noise.dat', format='ascii.basic')
    
    snr_bins = []
    for i in range(40):
            n_rec = len(match_r[(match_r/match_n>=2.75+0.5*i)*(match_r/match_n<3.25+0.5*i)])
            n_fake = len(fake_r[(fake_r/fake_n>=2.75+0.5*i)*(fake_r/fake_n<3.25+0.5*i)])
            n_sum = n_rec+n_fake
            if n_sum < 10: 
                if i >=5:
                    frac = 1.
                else: 
                    frac = 0.
            else: 
                frac = n_rec/n_sum
            snr_bins.append((0.5*i+3,frac,n_sum))
    snr, frac, num = zip(*snr_bins)

    t = Table([snr, frac, num], names=('snr', 'fidelity', 'number'))
    t.write(path[8:]+'fidelity_snr.dat', format='ascii.basic')
    print('N_spu='+str(np.mean(n_spu))+'+-'+str(np.std(n_spu)))
    return 0


if __name__ == '__main__': 
    err_max = 3.0 #float(input('The max of error: '))
    err_min = 0.5 #float(input('The min of error: '))

    flux_max = 26. #float(input('The max of flux: '))
    flux_min = 1. #float(input('The min of flux: '))

    path = '../mock_850/' #input('The name of the directory (e.g. mock_850): ') + '/'
    os.system('mkdir '+path[8:-1])
    n_trials = 10000 #int(input('The number of mock maps: '))
    #inik = int(input('The initial kernel size: '))
    #bigk = int(input('The number of different big kernels: '))

    #May have to change this based on how errors/fluxes distributed
    #same as correct_sources.py

    flux_step = (flux_max-flux_min)/25
    error_step = (err_max-err_min)/25
    
    flux_bins = np.linspace(flux_min, flux_max - flux_step, 25)
    error_bins = np.linspace(err_min, err_max - error_step, 25)
    
    '''max arcsecond separation between input and extracted
    source to qualify as a match'''
    thresh = 6. #int(input('The researching radius: '))
    analysis(n_trials)
    