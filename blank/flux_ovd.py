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

# calculate the (catalog & simulation) counts as function of flux in both inner and outer region, the flux bin is for every 2mJy
fsimin = []
fsimout = []
fsim = []
flux_bin = np.linspace(2,12,6)
flux = np.zeros(5)

# the mean value for each bin
for i in range(5):
    flux[i] = 0.5*flux_bin[i]+0.5*flux_bin[i+1]

# The number counts for each mock
for i in range(ntrials):
    fnc_sim = np.zeros(5)
    fnc_simin = np.zeros(5)
    fnc_simout = np.zeros(5)
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        
        for k in range(5):
            if flux_bin[k]<=f_sim[j]<flux_bin[k+1]:
                fnc_sim[k] = fnc_sim[k]+1
                
        if c_sim[j].separation(rg).arcmin<=4:
            for k in range(5):
                if flux_bin[k]<=f_sim[j]<flux_bin[k+1]:
                    fnc_simin[k] = fnc_simin[k]+1
        else:
            for k in range(5):
                if flux_bin[k]<=f_sim[j]<flux_bin[k+1]:
                    fnc_simout[k] = fnc_simout[k]+1
    fsim.append(fnc_sim)                
    fsimin.append(fnc_simin)
    fsimout.append(fnc_simout)


# The catalog number counts separate the inner and outer region
fn = np.zeros(5)
fnin = np.zeros(5)
fnout = np.zeros(5)
for i in range(len(c_my)):
    for j in range(5):
        if flux_bin[j]<=f_my[i]<flux_bin[j+1]:
            fn[j] = fn[j]+1

    if c_my[i].separation(rg).arcmin<=4:
        for j in range(5):
            if flux_bin[j]<=f_my[i]<flux_bin[j+1]:
                fnin[j] = fnin[j]+1
    else:
        for j in range(5):
            if flux_bin[j]<=f_my[i]<flux_bin[j+1]:
                fnout[j] = fnout[j]+1
                
# due to the low number (0) problem of the bright end, use the average of 20 maps
fsim = np.array(fsim)
fsimin = np.array(fsimin)
fsimout = np.array(fsimout)

np.save('flux/all.npy', np.array(fsim))
np.save('flux/in.npy', np.array(fsimin))
np.save('flux/out.npy', np.array(fsimout))

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
ov = []
ov_in = []
ov_out = []
ov_med = []
ov_in_med = []
ov_out_med = []
ov_84 = []
ov_in84 = []
ov_out84 = []
ov_16 = []
ov_in16 = []
ov_out16 = []

for i in range(int(ntrials/20)):
    fsimm = np.nanmean(fsim[i*20:(i+1)*20],axis=0)
    fsiminm = np.nanmean(fsimin[i*20:(i+1)*20],axis=0) # simulation inner region mean value
    fsimoutm = np.nanmean(fsimout[i*20:(i+1)*20],axis=0)
    #fin.append(fsiminm)
    #fout.append(fsimoutm)
    ov.append((fn-fsimm)/fsimm)
    ov_in.append((fnin-fsiminm)/fsiminm)
    ov_out.append((fnout-fsimoutm)/fsimoutm)

ov = np.array(ov)    
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
    ov_med.append(np.median(ov[:,i]))
    ov_in_med.append(np.median(ov_in[:,i]))
    ov_out_med.append(np.median(ov_out[:,i]))
    ov_84.append(np.nanpercentile(ov[:,i],84.1))
    ov_in84.append(np.nanpercentile(ov_in[:,i],84.1))
    ov_out84.append(np.nanpercentile(ov_out[:,i],84.1))
    ov_16.append(np.nanpercentile(ov[:,i],15.9))
    ov_in16.append(np.nanpercentile(ov_in[:,i],15.9))
    ov_out16.append(np.nanpercentile(ov_out[:,i],15.9))

ov_med = np.array(ov_med)
ov_in_med = np.array(ov_in_med)
ov_out_med = np.array(ov_out_med)
ov_84 = np.array(ov_84)
ov_in84 = np.array(ov_in84)
ov_out84 = np.array(ov_out84)
ov_16 = np.array(ov_16)
ov_in16 = np.array(ov_in16)
ov_out16 = np.array(ov_out16)    

# make correction due to 20 maps
#ov_84 = np.sqrt(20)*(ov_84-ov_med)/(ov_med+1)**2/fn + ov_med
#ov_in84 = np.sqrt(20)*(ov_in84-ov_in_med)/(ov_in_med+1)**2/fnin + ov_in_med
#ov_out84 = np.sqrt(20)*(ov_out84-ov_out_med)/(ov_out_med+1)**2/fnout + ov_out_med
#ov_16 = np.sqrt(20)*(ov_16-ov_med)/(ov_med+1)**2/fn + ov_med
#ov_in16 = np.sqrt(20)*(ov_in16-ov_in_med)/(ov_in_med+1)**2/fnin + ov_in_med
#ov_out16 = np.sqrt(20)*(ov_out16-ov_out_med)/(ov_out_med+1)**2/fnout + ov_out_med


# plot as a function of flux
plt.plot(flux,ov_med, color='tab:blue',label=r"Overdensity")
plt.plot(flux,ov_in_med, color='tab:green',label=r"Overdensity($\leq$4')")
plt.plot(flux,ov_out_med,color= 'tab:red',label="Overdensity(>4')" )
plt.scatter(flux,ov_med, color='tab:blue')
plt.scatter(flux,ov_in_med, color='tab:green')#,label="Overdensity(<4')")
plt.scatter(flux,ov_out_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.fill_between(flux,ov_16,ov_84,color='tab:blue',alpha=0.1 )
plt.fill_between(flux,ov_in16,ov_in84,color='tab:green',alpha=0.1 )
plt.fill_between(flux,ov_out16,ov_out84,color='tab:red',alpha=0.1 )

plt.ylim(-0.6,7.8)
plt.legend()
plt.xlabel('Flux (mJy)')
plt.ylabel('Overdensity')
plt.savefig('plot/ovd_f.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_f.eps',bbox_inches='tight')
plt.close()

# bar 
plt.bar(flux,ov_med,1.6,color='tab:blue',alpha = 0.6, label=r"Overdensity") # bar demonstration
plt.bar(flux-0.4,ov_in_med,0.8,color= 'tab:green',alpha = 0.6, label=r"Overdensity($\leq$4')" )
plt.bar(flux+0.4,ov_out_med,0.8,color= 'tab:red',alpha = 0.6, label="Overdensity(>4')" )
#plt.tick_params(labelsize = 15)
plt.legend()
plt.xlabel('Flux (mJy)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_f_bar.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_f_bar.eps',bbox_inches='tight')
plt.close()   


