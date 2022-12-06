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

# The real counts in different anuulus
ncr = np.zeros(9)
ncr_bright = np.zeros(9)
ncr_faint = np.zeros(9)
flux = []
flux_bright = []
flux_faint = []
for i in range(9):
    r_bin = []
    rb_bin = []
    rf_bin = []
    for j in range(len(c_my)):
        if i<=c_my[j].separation(rg).arcmin<i+1:
            ncr[i] = ncr[i]+1
            r_bin.append(f_my[j])
            if f_my[j]>=5: # from Garcia+2020, only the bright SMGs>5mJy would trace the massive structures at z~2
                ncr_bright[i] = ncr_bright[i]+1
                rb_bin.append(f_my[j])
            else:
                ncr_faint[i] = ncr_faint[i]+1
                rf_bin.append(f_my[j])

    flux.append(r_bin)
    flux_bright.append(rb_bin)
    flux_faint.append(rf_bin)


# The simulated counts, recall this would introduce a problem due to overdensity, need to use the average 
n_sim = []
nb_sim = []
nf_sim = []
for i in range(ntrials):
    nc_sim = np.zeros(9)
    ncb_sim = np.zeros(9) # bright sources >5mJy
    ncf_sim = np.zeros(9) # faint sources <5mJy
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(9):
            if k<=c_sim[j].separation(rg).arcmin<k+1:
                nc_sim[k] = nc_sim[k]+1
                if f_sim[j]>=5:
                    ncb_sim[k] = ncb_sim[k]+1
                else:
                    ncf_sim[k] = ncf_sim[k]+1
    n_sim.append(nc_sim)
    nb_sim.append(ncb_sim)
    nf_sim.append(ncf_sim)
    #ovd.append((ncr-nc_sim)/nc_sim) # unfinished

n_sim = np.array(n_sim)
nb_sim = np.array(nb_sim) 
nf_sim = np.array(nf_sim) 
np.save('radial/all.npy', np.array(n_sim))
np.save('radial/bright.npy', np.array(nb_sim))
np.save('radial/faint.npy', np.array(nf_sim))
# calculate the mean for every 20 maps
ov = []
ov_faint = []
ov_bright = []
ov_med = []
ov_faint_med = []
ov_bright_med = []
ov_84 = []
ov_faint_84 = []
ov_bright_84 = []
ov_16 = []
ov_faint_16 = []
ov_bright_16 = []
for i in range(int(ntrials/20)):
    n_sim_mean = np.nanmean(n_sim[i*20:(i+1)*20],axis=0) # each annulus mean value
    nb_sim_mean = np.nanmean(nb_sim[i*20:(i+1)*20],axis=0)
    nf_sim_mean = np.nanmean(nf_sim[i*20:(i+1)*20],axis=0)
    ov.append((ncr-n_sim_mean)/n_sim_mean)
    ov_bright.append((ncr_bright-nb_sim_mean)/nb_sim_mean)
    ov_faint.append((ncr_faint-nf_sim_mean)/nf_sim_mean)

ov = np.array(ov)
ov_bright = np.array(ov_bright)
ov_faint = np.array(ov_faint)
for i in range(9):
    ov_med.append(np.median(ov[:,i]))
    ov_bright_med.append(np.median(ov_bright[:,i]))
    ov_faint_med.append(np.median(ov_faint[:,i]))
    ov_84.append(np.nanpercentile(ov[:,i],84.1))
    ov_bright_84.append(np.nanpercentile(ov_bright[:,i],84.1))
    ov_faint_84.append(np.nanpercentile(ov_faint[:,i],84.1))
    ov_16.append(np.nanpercentile(ov[:,i],15.9))
    ov_bright_16.append(np.nanpercentile(ov_bright[:,i],15.9))
    ov_faint_16.append(np.nanpercentile(ov_faint[:,i],15.9))
    
ov_med = np.array(ov_med)
ov_faint_med = np.array(ov_faint_med)
ov_bright_med = np.array(ov_bright_med)
ov_84 = np.array(ov_84)
ov_faint_84 = np.array(ov_faint_84)
ov_bright_84 = np.array(ov_bright_84)
ov_16 = np.array(ov_16)
ov_faint_16 = np.array(ov_faint_16)
ov_bright_16 = np.array(ov_bright_16)

# correct for the 20 average
#ov_84 = np.sqrt(20)*(ov_84-ov_med)/(ov_med+1)**2/ncr + ov_med
#ov_bright_84 = np.sqrt(20)*(ov_bright_84-ov_bright_med)/(ov_bright_med+1)**2/ncr_bright + ov_bright_med
#ov_faint_84 = np.sqrt(20)*(ov_faint_84-ov_faint_med)/(ov_faint_med+1)**2/ncr_faint + ov_faint_med
#ov_16 = np.sqrt(20)*(ov_16-ov_med)/(ov_med+1)**2/ncr + ov_med
#ov_bright_16 = np.sqrt(20)*(ov_bright_16-ov_bright_med)/(ov_bright_med+1)**2/ncr_bright + ov_bright_med
#ov_faint_16 = np.sqrt(20)*(ov_faint_16-ov_faint_med)/(ov_faint_med+1)**2/ncr_faint + ov_faint_med

# plot as a function of radius
r = np.linspace(0.5,8.5,9)
plt.plot(r,ov_med, color='tab:blue',label="Overdensity")
plt.plot(r,ov_bright_med,color= 'tab:red',label=r"Overdensity($\geq$5mJy)" )
plt.plot(r,ov_faint_med,color= 'tab:green',label="Overdensity(<5mJy)" )
plt.scatter(r,ov_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(r,ov_bright_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.scatter(r,ov_faint_med,color= 'tab:green')#,label="Overdensity(>4')" )
plt.fill_between(r,ov_16,ov_84,color='tab:blue',alpha=0.1 )
plt.fill_between(r,ov_bright_16,ov_bright_84,color='tab:red',alpha=0.1 )
plt.fill_between(r,ov_faint_16,ov_faint_84,color='tab:green',alpha=0.1 )

plt.legend()
plt.xlabel('Distance (arcmin)')
plt.ylabel('Overdensity')
plt.savefig('plot/ovd_r.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_r.eps',bbox_inches='tight')
plt.close()

# bar 
plt.bar(r,ov_med,0.8,color='tab:blue',alpha = 0.6, label="Overdensity") # bar demonstration
plt.bar(r+0.2,ov_bright_med,0.4,color= 'tab:red',alpha = 0.6, label=r"Overdensity($\geq$5mJy)" )
plt.bar(r-0.2,ov_faint_med,0.4,color= 'tab:green',alpha = 0.6, label=r"Overdensity(<5mJy)" )
#plt.tick_params(labelsize = 15)
plt.legend()
plt.xlabel('Distance (arcmin)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_r_bar.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_r_bar.eps',bbox_inches='tight')
plt.close()   




# overall overdensity calculation
n_spu = []
n_rec = []
nb_rec = []
nf_rec = []
ovo_whole = [] # the whole overdensity
ovo_bright = [] # the overdensity only for bright sources 
ovo_faint = []
for i in range(ntrials):
    fr = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')[0]
    n_rec.append(len(fr))
    nb_rec.append(len(fr[fr>=5]))
    nf_rec.append(len(fr[fr<5]))
    ovo_whole.append((len(f_my)-len(fr))/len(fr))
    ovo_bright.append((len(f_my[f_my>=5])-len(fr[fr>=5]))/len(fr[fr>=5]))
    ovo_faint.append((len(f_my[f_my<5])-len(fr[fr<5]))/len(fr[fr<5]))
    try:
        fs = read_spu('mock_850/spu_map'+str(i)+'.dat')[0]
        n_spu.append(len(fs))
    except:
        n_spu.append(0)
        
#mean_ovd = (len(f_my)-np.nanmean(n_rec))/np.nanmean(n_rec)
np.save('all/rec.npy', np.array(n_rec))
np.save('all/spu.npy', np.array(n_spu))
np.save('all/bright.npy', np.array(nb_rec))
np.save('all/faint.npy', np.array(nf_rec))
print('The number of spurious sources and recovered sources (also bright sources >=5mJy only and faint sources < 5mJy) and their corresponding standard deviations in the simulation: \n',np.nanmean(n_spu),np.nanstd(n_spu),np.nanmean(n_rec),np.nanstd(n_rec),np.nanmean(nb_rec),np.nanstd(nb_rec),np.nanmean(nf_rec),np.nanstd(nf_rec))
print('The overall overdensity: ', np.nanmean(ovo_whole),'±',np.nanstd(ovo_whole))
print('The overdensity for s >= 5 mJy: ', np.nanmean(ov_bright),'±',np.nanstd(ov_bright))
print('The overdensity for s < 5 mJy: ', np.nanmean(ov_faint),'±',np.nanstd(ov_faint))

