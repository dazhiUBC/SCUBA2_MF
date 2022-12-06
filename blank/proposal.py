import pandas as pd
from MF import *
from astropy.coordinates import SkyCoord

# this is just for the proposal, don't use it
ntrials = 500 # number of mock maps

def read_ca(fname):
    cata = pd.read_csv(fname,delimiter= ' ' )
    f = cata['flux'] # flux
    e = cata['err'] # err
    c = SkyCoord(ra=cata['ra'],dec=cata['dec'],unit=(u.hourangle, u.deg)) # coordinate
    return f,e,c

def read_sim(fname):
    cata = pd.read_csv(fname,delimiter= ' ' )
    f = cata['flux']
    c = SkyCoord(ra=cata['ra'],dec=cata['dec'],unit=(u.hourangle, u.deg))
    return f,c

def read_spu(fname):
    cata = pd.read_csv(fname,delimiter= ' ' )
    f = cata['spu_flux']
    c = SkyCoord(ra=cata['ra'],dec=cata['dec'],unit=(u.hourangle, u.deg))
    return f,c

# read the actual catalog
f_my, e_my, c_my  = read_ca('../sources_4C23_850_cal_crop_MF.dat')

# define the position of the HzRGs
rg = SkyCoord(ra=316.811766,dec=23.529172, unit= u.deg)

# The real counts in different anuulus
ncr = np.zeros(9)
ncr_bright = np.zeros(9)
flux = []
flux_bright = []
for i in range(9):
    r_bin = []
    rb_bin = []
    for j in range(len(c_my)):
        if i<=c_my[j].separation(rg).arcmin<i+1:
            ncr[i] = ncr[i]+1
            r_bin.append(f_my[j])
            if f_my[j]>5: # from Garcia+2020, only the bright SMGs>5mJy would trace the massive structures at z~2
                ncr_bright[i] = ncr_bright[i]+1
                rb_bin.append(f_my[j])
    flux.append(r_bin)
    flux_bright.append(rb_bin)


# The simulated counts, recall this would introduce a problem due to overdensity, need to use the average 
n_sim = []
nb_sim = []
for i in range(ntrials):
    nc_sim = np.zeros(9)
    ncb_sim = np.zeros(9) # bright sources >5mJy
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(9):
            if k<=c_sim[j].separation(rg).arcmin<k+1:
                nc_sim[k] = nc_sim[k]+1
                if f_sim[j]>5:
                    ncb_sim[k] = ncb_sim[k]+1
    n_sim.append(nc_sim)
    nb_sim.append(ncb_sim)
    #ovd.append((ncr-nc_sim)/nc_sim) # unfinished

n_sim = np.array(n_sim)
nb_sim = np.array(nb_sim)    

# calculate the mean for every 20 maps
ov = []
ov_bright = []
ov_med = []
ov_bright_med = []
ov_68 = []
ov_bright_68 = []
ov_32 = []
ov_bright_32 = []
for i in range(int(ntrials/20)):
    n_sim_mean = np.nanmean(n_sim[i*20:(i+1)*20],axis=0) # each annulus mean value
    nb_sim_mean = np.nanmean(nb_sim[i*20:(i+1)*20],axis=0)
    ov.append((ncr-n_sim_mean)/n_sim_mean)
    ov_bright.append((ncr_bright-nb_sim_mean)/nb_sim_mean)

ov = np.array(ov)
ov_bright = np.array(ov_bright)

for i in range(9):
    ov_med.append(np.median(ov[:,i]))
    ov_bright_med.append(np.median(ov_bright[:,i]))
    ov_68.append(np.nanpercentile(ov[:,i],68))
    ov_bright_68.append(np.nanpercentile(ov_bright[:,i],68))
    ov_32.append(np.nanpercentile(ov[:,i],32))
    ov_bright_32.append(np.nanpercentile(ov_bright[:,i],32))

# plot as a function of radius
r = np.linspace(0.5,8.5,9)
plt.plot(r,ov_med, color='tab:blue',label="Overdensity")
plt.plot(r,ov_bright_med,color= 'tab:red',label="Overdensity(>5mJy)" )
plt.scatter(r,ov_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(r,ov_bright_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.fill_between(r,ov_32,ov_68,color='tab:blue',alpha=0.1 )
plt.fill_between(r,ov_bright_32,ov_bright_68,color='tab:red',alpha=0.1 )

plt.legend()
plt.xlabel('Distance (arcmin)')
plt.ylabel('Overdensity')
plt.savefig('proposal/ovd_r.pdf')
plt.savefig('proposal/ovd_r.eps')
plt.close()



# overall overdensity calculation
n_spu = []
n_rec = []
nb_rec = []
ov_whole = [] # the whole overdensity
ov_bright = [] # the overdensity only for bright sources 
for i in range(ntrials):
    fr = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')[0]
    n_rec.append(len(fr))
    nb_rec.append(len(fr[fr>5]))
    ov_whole.append((len(f_my)-len(fr))/len(fr))
    ov_bright.append((len(f_my[f_my>5])-len(fr[fr>5]))/len(fr[fr>5]))
    try:
        fs = read_spu('mock_850/spu_map'+str(i)+'.dat')[0]
        n_spu.append(len(fs))
    except:
        n_spu.append(0)
        
#mean_ovd = (len(f_my)-np.nanmean(n_rec))/np.nanmean(n_rec)

print('The number of spurious sources and recovered sources (also bright sources >5mJy only) and their corresponding standard deviations in the simulation: \n',np.nanmean(n_spu),np.nanstd(n_spu),np.nanmean(n_rec),np.nanstd(n_rec),np.nanmean(nb_rec),np.nanstd(nb_rec))
print('The overall overdensity: ', np.nanmean(ov_whole),'±',np.nanstd(ov_whole))
print('The overdensity for s > 5 mJy: ', np.nanmean(ov_bright),'±',np.nanstd(ov_bright))




# calculate the (catalog & simulation) counts as function of flux in both inner and outer region, the flux bin is for every 2mJy
fsimin = []
fsimout = []
flux_bin = np.linspace(2,12,6)
flux = np.zeros(5)

# the mean value for each bin
for i in range(5):
    flux[i] = 0.5*flux_bin[i]+0.5*flux_bin[i+1]

# The number counts for each mock
for i in range(ntrials):
    fnc_sim = np.zeros(5)
    fnc_simout = np.zeros(5)
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        if c_sim[j].separation(rg).arcmin<=4:
            for k in range(5):
                if flux_bin[k]<=f_sim[j]<flux_bin[k+1]:
                    fnc_sim[k] = fnc_sim[k]+1
        else:
            for k in range(5):
                if flux_bin[k]<=f_sim[j]<flux_bin[k+1]:
                    fnc_simout[k] = fnc_simout[k]+1
                    
    fsimin.append(fnc_sim)
    fsimout.append(fnc_simout)


# The catalog number counts separate the inner and outer region
fn = np.zeros(5)
fnout = np.zeros(5)
for i in range(len(c_my)):
    if c_my[i].separation(rg).arcmin<=4:
        for j in range(5):
            if flux_bin[j]<=f_my[i]<flux_bin[j+1]:
                fn[j] = fn[j]+1
    else:
        for j in range(5):
            if flux_bin[j]<=f_my[i]<flux_bin[j+1]:
                fnout[j] = fnout[j]+1
                
# due to the low number (0) problem of the bright end, use the average of 20 maps

fsimin = np.array(fsimin)
fsimout = np.array(fsimout)

#f_min = []
#f_sin = [] it is not gaussian distribution, therefore standard deviation not accurate
#f_mout = []
#f_sout = []
#fin = []
#fout = []
#fin_min = [] # 68 and 32
#fin_max = []
#fout_min = []
#fout_max = []
ov_in = []
ov_out = []
ov_in_med = []
ov_out_med = []
ov_in68 = []
ov_out68 = []
ov_in32 = []
ov_out32 = []

for i in range(int(ntrials/20)):
    fsiminm = np.nanmean(fsimin[i*20:(i+1)*20],axis=0) # simulation inner region mean value
    fsimoutm = np.nanmean(fsimout[i*20:(i+1)*20],axis=0)
    #fin.append(fsiminm)
    #fout.append(fsimoutm)
    ov_in.append((fn-fsiminm)/fsiminm)
    ov_out.append((fnout-fsimoutm)/fsimoutm)

ov_in = np.array(ov_in)
ov_out = np.array(ov_out)

# statistic
for i in range(5):
    #f_min.append(np.median(fin[:,i]))
    #f_sin.append(np.nanstd(fin[:,i]))
    #f_mout.append(np.median(fout[:,i]))
    #f_sout.append(np.nanstd(fout[:,i]))
    #fin_min.append(np.nanpercentile(fin[:,i],32))
    #fin_max.append(np.nanpercentile(fin[:,i],68))
    #fout_min.append(np.nanpercentile(fout[:,i],32))
    #fout_max.append(np.nanpercentile(fout[:,i],68))
    ov_in_med.append(np.median(ov_in[:,i]))
    ov_out_med.append(np.median(ov_out[:,i]))
    ov_in68.append(np.nanpercentile(ov_in[:,i],68))
    ov_out68.append(np.nanpercentile(ov_out[:,i],68))
    ov_in32.append(np.nanpercentile(ov_in[:,i],32))
    ov_out32.append(np.nanpercentile(ov_out[:,i],32))
    
# plot as a function of flux
plt.plot(flux,ov_in_med, color='tab:blue',label="Overdensity(<4')")
plt.plot(flux,ov_out_med,color= 'tab:red',label="Overdensity(>4')" )
plt.scatter(flux,ov_in_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(flux,ov_out_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.fill_between(flux,ov_in32,ov_in68,color='tab:blue',alpha=0.1 )

plt.fill_between(flux,ov_out32,ov_out68,color='tab:red',alpha=0.1 )

#plt.ylim(-0.6,6.5)
plt.legend()
plt.xlabel('Flux (mJy)')
plt.ylabel('Overdensity')
plt.savefig('proposal/ovd_f.pdf')
plt.savefig('proposal/ovd_f.eps')
plt.close()



# finally, plot the overdensity as a function of angular position (for every 30 degree)
# of course, it can be combined with the radial and flux analysis and reduce the time.
# This script is based on multiple scripts that I used in the past and I am too lazy to modify it.
# Hopefully it won't take too much time... (first rule: if it works, don't touch it)
# as other suggested, it would be plotted as historgram
nca = np.zeros(12) # number count for different angular position
nca_bright = np.zeros(12)
for i in range(12):
    for j in range(len(c_my)):
        if 30*i<=c_my[j].position_angle(rg).degree<30*(i+1):
            nca[i] = nca[i]+1
            if f_my[j]>5: 
                nca_bright[i] = nca_bright[i]+1

na_sim = []
nab_sim = []
for i in range(ntrials):
    nca_sim = np.zeros(12)
    ncab_sim = np.zeros(12) # bright sources > 5 mJy
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(12):
            if 30*k<=c_sim[j].position_angle(rg).degree<30*(k+1):
                nca_sim[k] = nca_sim[k]+1
                if f_sim[j]>5:
                    ncab_sim[k] = ncab_sim[k]+1
    na_sim.append(nca_sim)
    nab_sim.append(ncab_sim)

na_sim = np.array(na_sim)
nab_sim = np.array(nab_sim)    

# calculate the mean for every 20 maps
ova = []
ova_bright = []
ova_med = []
ova_bright_med = []
ova_68 = []
ova_bright_68 = []
ova_32 = []
ova_bright_32 = []
for i in range(int(ntrials/20)):
    na_sim_mean = np.nanmean(na_sim[i*20:(i+1)*20],axis=0) # each annulus mean value
    nab_sim_mean = np.nanmean(nab_sim[i*20:(i+1)*20],axis=0)
    ova.append((nca-na_sim_mean)/na_sim_mean)
    ova_bright.append((nca_bright-nab_sim_mean)/nab_sim_mean)

ova = np.array(ova)
ova_bright = np.array(ova_bright)

for i in range(12):
    ova_med.append(np.median(ova[:,i]))
    ova_bright_med.append(np.median(ova_bright[:,i]))
    ova_68.append(np.nanpercentile(ova[:,i],68))
    ova_bright_68.append(np.nanpercentile(ova_bright[:,i],68))
    ova_32.append(np.nanpercentile(ova[:,i],32))
    ova_bright_32.append(np.nanpercentile(ova_bright[:,i],32))

# plot as a function of radius
a = np.linspace(15,345,12)
plt.plot(a,ova_med, color='tab:blue',label="Overdensity")
plt.plot(a,ova_bright_med,color= 'tab:red',label="Overdensity(>5mJy)" )
plt.scatter(a,ova_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(a,ova_bright_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.fill_between(a,ova_32,ova_68,color='tab:blue',alpha=0.1 )
plt.fill_between(a,ova_bright_32,ova_bright_68,color='tab:red',alpha=0.1 )

plt.legend()
plt.xlabel('Position angle (degree)')
plt.ylabel('Overdensity')
plt.savefig('proposal/ovd_a.pdf')
plt.savefig('proposal/ovd_a.eps') # traditional visualization
plt.close()

plt.bar(a-7,ova_med,14,color='tab:blue',label="Overdensity") # bar demonstration
plt.bar(a+7,ova_bright_med,14,color= 'tab:red',label="Overdensity(>5mJy)" )

plt.legend()
plt.xlabel('Position angle (degree)')
plt.ylabel('Overdensity')
plt.savefig('proposal/ovd_a_bar.pdf')
plt.savefig('proposal/ovd_a_bar.eps')
plt.close()                