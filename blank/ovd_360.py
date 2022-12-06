import pandas as pd
from MF import *
from astropy.coordinates import SkyCoord

ntrials = 10000 # number of mock maps

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

# finally, plot the overdensity as a function of angular position (for every 30 degree)
# of course, it can be combined with the radial and flux analysis and reduce the time.
# This script is based on multiple scripts that I used in the past and I am too lazy to modify it.
# Hopefully it won't take too much time... (first rule: if it works, don't touch it)
# as other suggested, it would be plotted as historgram
nca = np.zeros(12) # number count for different angular position
nca_bright = np.zeros(12)
nca_faint = np.zeros(12)
for i in range(12):
    for j in range(len(c_my)):
        if 30*i<=c_my[j].position_angle(rg).degree<30*(i+1):
            nca[i] = nca[i]+1
            if f_my[j]>=5: 
                nca_bright[i] = nca_bright[i]+1
            else:
                nca_faint[i] = nca_faint[i]+1

na_sim = []
nab_sim = []
naf_sim = []
for i in range(ntrials):
    nca_sim = np.zeros(12)
    ncab_sim = np.zeros(12) # bright sources > 5 mJy
    ncaf_sim = np.zeros(12)
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(12):
            if 30*k<=c_sim[j].position_angle(rg).degree<30*(k+1):
                nca_sim[k] = nca_sim[k]+1
                if f_sim[j]>=5:
                    ncab_sim[k] = ncab_sim[k]+1
                else:
                    ncaf_sim[k] = ncaf_sim[k]+1
    na_sim.append(nca_sim)
    nab_sim.append(ncab_sim)
    naf_sim.append(ncaf_sim)

na_sim = np.array(na_sim)
nab_sim = np.array(nab_sim)    
naf_sim = np.array(naf_sim) 

np.save('angular360/all.npy', np.array(na_sim))
np.save('angular360/bright.npy', np.array(nab_sim))
np.save('angular360/faint.npy', np.array(naf_sim))

# calculate the mean for every 20 maps
ova = []
ova_bright = []
ova_med = []
ova_bright_med = []
ova_84 = []
ova_bright_84 = []
ova_16 = []
ova_bright_16 = []
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
    ova_84.append(np.nanpercentile(ova[:,i],84.1))
    ova_bright_84.append(np.nanpercentile(ova_bright[:,i],84.1))
    ova_16.append(np.nanpercentile(ova[:,i],15.9))
    ova_bright_16.append(np.nanpercentile(ova_bright[:,i],15.9))

# plot as a function of angle
a = np.linspace(15,345,12)
plt.plot(a,ova_med, color='tab:blue',label="Overdensity")
plt.plot(a,ova_bright_med,color= 'tab:red',label="Overdensity(>5mJy)" )
plt.scatter(a,ova_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(a,ova_bright_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.fill_between(a,ova_16,ova_84,color='tab:blue',alpha=0.1 )
plt.fill_between(a,ova_bright_16,ova_bright_84,color='tab:red',alpha=0.1 )

plt.legend()
plt.xlabel('Position angle (degree)')
plt.ylabel('Overdensity')
plt.savefig('plot/360ovd_a.pdf',bbox_inches='tight')
plt.savefig('plot/360ovd_a.eps',bbox_inches='tight') # traditional visualization
plt.close()

plt.bar(a-7,ova_med,14,color='tab:blue',label="Overdensity") # bar demonstration
plt.bar(a+7,ova_bright_med,14,color= 'tab:red',label="Overdensity(>5mJy)" )

plt.legend()
plt.xlabel('Position angle (degree)')
plt.ylabel('Overdensity')
plt.savefig('plot/360ovd_a_bar.pdf',bbox_inches='tight')
plt.savefig('plot/360ovd_a_bar.eps',bbox_inches='tight')
plt.close()                