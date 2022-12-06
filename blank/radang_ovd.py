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
ncr_jet = np.zeros(9)
ncr_bright_jet = np.zeros(9)
ncr_faint_jet = np.zeros(9)

ncr_per = np.zeros(9)
ncr_bright_per = np.zeros(9)
ncr_faint_per = np.zeros(9)

#flux = []
#flux_bright = []
#flux_faint = []
for i in range(9):
    #r_bin = []
    #rb_bin = []
    #rf_bin = []
    for j in range(len(c_my)):
        if 0<=c_my[j].position_angle(rg).degree<90 or 180<=c_my[j].position_angle(rg).degree<270:
            if i<=c_my[j].separation(rg).arcmin<i+1:
                ncr_jet[i] = ncr_jet[i]+1
                #r_bin.append(f_my[j])
                if f_my[j]>=5: # from Garcia+2020, only the bright SMGs>5mJy would trace the massive structures at z~2
                    ncr_bright_jet[i] = ncr_bright_jet[i]+1
                    #rb_bin.append(f_my[j])
                else:
                    ncr_faint_jet[i] = ncr_faint_jet[i]+1
                    #rf_bin.append(f_my[j])
        else:
            if i<=c_my[j].separation(rg).arcmin<i+1:
                ncr_per[i] = ncr_per[i]+1
                if f_my[j]>=5: 
                    ncr_bright_per[i] = ncr_bright_per[i]+1
                else:
                    ncr_faint_per[i] = ncr_faint_per[i]+1


    #flux.append(r_bin)
    #flux_bright.append(rb_bin)
    #flux_faint.append(rf_bin)


# The simulated counts, recall this would introduce a problem due to overdensity, need to use the average 
n_sim_jet = []
nb_sim_jet = []
nf_sim_jet = []

n_sim_per = []
nb_sim_per = []
nf_sim_per = []

for i in range(ntrials):
    nc_sim_jet = np.zeros(9)
    ncb_sim_jet = np.zeros(9) # bright sources >5mJy
    ncf_sim_jet = np.zeros(9) # faint sources <5mJy
    
    nc_sim_per = np.zeros(9)
    ncb_sim_per = np.zeros(9) # bright sources >5mJy
    ncf_sim_per = np.zeros(9) # faint sources <5mJy

    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(9):
            if 0<=c_sim[j].position_angle(rg).degree<90 or 180<=c_sim[j].position_angle(rg).degree<270:
                if k<=c_sim[j].separation(rg).arcmin<k+1:
                    nc_sim_jet[k] = nc_sim_jet[k]+1
                    if f_sim[j]>=5:
                        ncb_sim_jet[k] = ncb_sim_jet[k]+1
                    else:
                        ncf_sim_jet[k] = ncf_sim_jet[k]+1
            else:
                if k<=c_sim[j].separation(rg).arcmin<k+1:
                    nc_sim_per[k] = nc_sim_per[k]+1
                    if f_sim[j]>=5:
                        ncb_sim_per[k] = ncb_sim_per[k]+1
                    else:
                        ncf_sim_per[k] = ncf_sim_per[k]+1

    n_sim_jet.append(nc_sim_jet)
    nb_sim_jet.append(ncb_sim_jet)
    nf_sim_jet.append(ncf_sim_jet)
    n_sim_per.append(nc_sim_per)
    nb_sim_per.append(ncb_sim_per)
    nf_sim_per.append(ncf_sim_per)
    #ovd.append((ncr-nc_sim)/nc_sim) # unfinished

n_sim_jet = np.array(n_sim_jet)
nb_sim_jet = np.array(nb_sim_jet) 
nf_sim_jet = np.array(nf_sim_jet) 

n_sim_per = np.array(n_sim_per)
nb_sim_per = np.array(nb_sim_per) 
nf_sim_per = np.array(nf_sim_per) 

np.save('angrad/all_per.npy', np.array(n_sim_per))
np.save('angrad/bright_per.npy', np.array(nb_sim_per))
np.save('angrad/faint_per.npy', np.array(nf_sim_per))

np.save('angrad/all_jet.npy', np.array(n_sim_jet))
np.save('angrad/bright_jet.npy', np.array(nb_sim_jet))
np.save('angrad/faint_jet.npy', np.array(nf_sim_jet))

# calculate the mean for every 20 maps
ov_jet = []
ov_faint_jet = []
ov_bright_jet = []
ov_med_jet = []
ov_faint_med_jet = []
ov_bright_med_jet = []
ov_84_jet = []
ov_faint_84_jet = []
ov_bright_84_jet = []
ov_16_jet = []
ov_faint_16_jet = []
ov_bright_16_jet = []

ov_per = []
ov_faint_per = []
ov_bright_per = []
ov_med_per = []
ov_faint_med_per = []
ov_bright_med_per = []
ov_84_per = []
ov_faint_84_per = []
ov_bright_84_per = []
ov_16_per = []
ov_faint_16_per = []
ov_bright_16_per = []


for i in range(int(ntrials/25)):
    n_sim_mean_jet = np.nanmean(n_sim_jet[i*25:(i+1)*25],axis=0) # each annulus mean value
    nb_sim_mean_jet = np.nanmean(nb_sim_jet[i*25:(i+1)*25],axis=0)
    nf_sim_mean_jet = np.nanmean(nf_sim_jet[i*25:(i+1)*25],axis=0)
    ov_jet.append((ncr_jet-n_sim_mean_jet)/n_sim_mean_jet)
    ov_bright_jet.append((ncr_bright_jet-nb_sim_mean_jet)/nb_sim_mean_jet)
    ov_faint_jet.append((ncr_faint_jet-nf_sim_mean_jet)/nf_sim_mean_jet)
    
    n_sim_mean_per = np.nanmean(n_sim_per[i*25:(i+1)*25],axis=0) # each annulus mean value
    nb_sim_mean_per = np.nanmean(nb_sim_per[i*25:(i+1)*25],axis=0)
    nf_sim_mean_per = np.nanmean(nf_sim_per[i*25:(i+1)*25],axis=0)
    ov_per.append((ncr_per-n_sim_mean_per)/n_sim_mean_per)
    ov_bright_per.append((ncr_bright_per-nb_sim_mean_per)/nb_sim_mean_per)
    ov_faint_per.append((ncr_faint_per-nf_sim_mean_per)/nf_sim_mean_per)


ov_jet = np.array(ov_jet)
ov_bright_jet = np.array(ov_bright_jet)
ov_faint_jet = np.array(ov_faint_jet)

ov_per = np.array(ov_per)
ov_bright_per = np.array(ov_bright_per)
ov_faint_per = np.array(ov_faint_per)

for i in range(9):
    ov_med_jet.append(np.median(ov_jet[:,i]))
    ov_bright_med_jet.append(np.median(ov_bright_jet[:,i]))
    ov_faint_med_jet.append(np.median(ov_faint_jet[:,i]))
    ov_84_jet.append(np.nanpercentile(ov_jet[:,i],84.1))
    ov_bright_84_jet.append(np.nanpercentile(ov_bright_jet[:,i],84.1))
    ov_faint_84_jet.append(np.nanpercentile(ov_faint_jet[:,i],84.1))
    ov_16_jet.append(np.nanpercentile(ov_jet[:,i],15.9))
    ov_bright_16_jet.append(np.nanpercentile(ov_bright_jet[:,i],15.9))
    ov_faint_16_jet.append(np.nanpercentile(ov_faint_jet[:,i],15.9))
    
    ov_med_per.append(np.median(ov_per[:,i]))
    ov_bright_med_per.append(np.median(ov_bright_per[:,i]))
    ov_faint_med_per.append(np.median(ov_faint_per[:,i]))
    ov_84_per.append(np.nanpercentile(ov_per[:,i],84.1))
    ov_bright_84_per.append(np.nanpercentile(ov_bright_per[:,i],84.1))
    ov_faint_84_per.append(np.nanpercentile(ov_faint_per[:,i],84.1))
    ov_16_per.append(np.nanpercentile(ov_per[:,i],15.9))
    ov_bright_16_per.append(np.nanpercentile(ov_bright_per[:,i],15.9))
    ov_faint_16_per.append(np.nanpercentile(ov_faint_per[:,i],15.9))


    
ov_med_jet = np.array(ov_med_jet)
ov_faint_med_jet = np.array(ov_faint_med_jet)
ov_bright_med_jet = np.array(ov_bright_med_jet)
ov_84_jet = np.array(ov_84_jet)
ov_faint_84_jet = np.array(ov_faint_84_jet)
ov_bright_84_jet = np.array(ov_bright_84_jet)
ov_16_jet = np.array(ov_16_jet)
ov_faint_16_jet = np.array(ov_faint_16_jet)
ov_bright_16_jet = np.array(ov_bright_16_jet)

ov_med_per = np.array(ov_med_per)
ov_faint_med_per = np.array(ov_faint_med_per)
ov_bright_med_per = np.array(ov_bright_med_per)
ov_84_per = np.array(ov_84_per)
ov_faint_84_per = np.array(ov_faint_84_per)
ov_bright_84_per = np.array(ov_bright_84_per)
ov_16_per = np.array(ov_16_per)
ov_faint_16_per = np.array(ov_faint_16_per)
ov_bright_16_per = np.array(ov_bright_16_per)

# correct for the 25 average
ov_84_jet = 5*(ov_84_jet-ov_med_jet) + ov_med_jet
ov_bright_84_jet = 5*(ov_bright_84_jet-ov_bright_med_jet) + ov_bright_med_jet
ov_faint_84_jet =5*(ov_faint_84_jet-ov_faint_med_jet) + ov_faint_med_jet
ov_16_jet = 5*(ov_16_jet-ov_med_jet) + ov_med_jet
ov_bright_16_jet = 5*(ov_bright_16_jet-ov_bright_med_jet) + ov_bright_med_jet
ov_faint_16_jet = 5*(ov_faint_16_jet-ov_faint_med_jet) + ov_faint_med_jet

ov_84_per = 5*(ov_84_per-ov_med_per) + ov_med_per
ov_bright_84_per = 5*(ov_bright_84_per-ov_bright_med_per) + ov_bright_med_per
ov_faint_84_per =5*(ov_faint_84_per-ov_faint_med_per) + ov_faint_med_per
ov_16_per = 5*(ov_16_per-ov_med_per) + ov_med_per
ov_bright_16_per = 5*(ov_bright_16_per-ov_bright_med_per) + ov_bright_med_per
ov_faint_16_per = 5*(ov_faint_16_per-ov_faint_med_per) + ov_faint_med_per

# plot as a function of radius
r = np.linspace(0.5,8.5,9)
plt.plot(r,ov_med_jet, color='tab:blue',label="Overdensity(jet)")
plt.plot(r,ov_bright_med_jet,color= 'tab:red',label=r"Overdensity($\geq$5mJy)(jet)" )
plt.plot(r,ov_faint_med_jet,color= 'tab:green',label="Overdensity(<5mJy)(jet)" )
plt.scatter(r,ov_med_jet, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(r,ov_bright_med_jet,color= 'tab:red')#,label="Overdensity(>4')" )
plt.scatter(r,ov_faint_med_jet,color= 'tab:green')#,label="Overdensity(>4')" )
plt.fill_between(r,ov_16_jet,ov_84_jet,color='tab:blue',alpha=0.1 )
plt.fill_between(r,ov_bright_16_jet,ov_bright_84_jet,color='tab:red',alpha=0.1 )
plt.fill_between(r,ov_faint_16_jet,ov_faint_84_jet,color='tab:green',alpha=0.1 )
plt.ylim(-1.1,4.9)
plt.legend()
plt.xlabel('Distance (arcmin)')
plt.ylabel('Overdensity')
plt.savefig('plot/ovd_r_jet.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_r_jet.eps',bbox_inches='tight')
plt.close()

plt.plot(r,ov_med_per, color='tab:blue',label="Overdensity(per)")
plt.plot(r,ov_bright_med_per,color= 'tab:red',label=r"Overdensity($\geq$5mJy)(per)" )
plt.plot(r,ov_faint_med_per,color= 'tab:green',label="Overdensity(<5mJy)(per)" )
plt.scatter(r,ov_med_per, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(r,ov_bright_med_per,color= 'tab:red')#,label="Overdensity(>4')" )
plt.scatter(r,ov_faint_med_per,color= 'tab:green')#,label="Overdensity(>4')" )
plt.fill_between(r,ov_16_per,ov_84_per,color='tab:blue',alpha=0.1 )
plt.fill_between(r,ov_bright_16_per,ov_bright_84_per,color='tab:red',alpha=0.1 )
plt.fill_between(r,ov_faint_16_per,ov_faint_84_per,color='tab:green',alpha=0.1 )
plt.ylim(-1.1,4.9)
plt.legend()
plt.xlabel('Distance (arcmin)')
plt.ylabel('Overdensity')
plt.savefig('plot/ovd_r_per.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_r_per.eps',bbox_inches='tight')
plt.close()


'''
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

'''