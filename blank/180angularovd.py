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

nca = np.zeros(12) # number count for different angular position
nca_bright = np.zeros(12)
nca_faint = np.zeros(12)
for i in range(12):
    for j in range(len(c_my)):
        if 15*i<=rg.position_angle(c_my[j]).degree<15*(i+1) or 180+15*i<=rg.position_angle(c_my[j]).degree<180+15*(i+1):
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
            if 15*k<=rg.position_angle(c_sim[j]).degree<15*(k+1) or 180+15*k<=rg.position_angle(c_sim[j]).degree<180+15*(k+1):
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
np.save('angular/all.npy', np.array(na_sim))
np.save('angular/bright.npy', np.array(nab_sim))
np.save('angular/faint.npy', np.array(naf_sim))


# calculate the mean for every 20 maps
ova = []
ova_bright = []
ova_faint = []
ova_med = []
ova_bright_med = []
ova_faint_med = []
ova_84 = []
ova_bright_84 = []
ova_faint_84 = []
ova_16 = []
ova_bright_16 = []
ova_faint_16 = []
for i in range(int(ntrials/4)):
    na_sim_mean = np.nanmean(na_sim[i*4:(i+1)*4],axis=0) # each annulus mean value
    nab_sim_mean = np.nanmean(nab_sim[i*4:(i+1)*4],axis=0)
    naf_sim_mean = np.nanmean(naf_sim[i*4:(i+1)*4],axis=0)
    ova.append((nca-na_sim_mean)/na_sim_mean)
    ova_bright.append((nca_bright-nab_sim_mean)/nab_sim_mean)
    ova_faint.append((nca_faint-naf_sim_mean)/naf_sim_mean)

ova = np.array(ova)
ova_bright = np.array(ova_bright)
ova_faint = np.array(ova_faint)

for i in range(12):
    ova_med.append(np.median(ova[:,i]))
    ova_bright_med.append(np.median(ova_bright[:,i]))
    ova_faint_med.append(np.median(ova_faint[:,i]))
    ova_84.append(np.nanpercentile(ova[:,i],84.1))
    ova_bright_84.append(np.nanpercentile(ova_bright[:,i],84.1))
    ova_faint_84.append(np.nanpercentile(ova_faint[:,i],84.1))
    ova_16.append(np.nanpercentile(ova[:,i],15.9))
    ova_bright_16.append(np.nanpercentile(ova_bright[:,i],15.9))
    ova_faint_16.append(np.nanpercentile(ova_faint[:,i],15.9))
ova_med = np.array(ova_med)
ova_bright_med = np.array(ova_bright_med)
ova_faint_med = np.array(ova_faint_med)
ova_84 = np.array(ova_84)
ova_bright_84 = np.array(ova_bright_84)
ova_faint_84 = np.array(ova_faint_84)
ova_16 = np.array(ova_16)
ova_bright_16 = np.array(ova_bright_16)
ova_faint_16 = np.array(ova_faint_16)
ova_84 = 3*(ova_84-ova_med) + ova_med
ova_bright_84 = 3*(ova_bright_84-ova_bright_med)+ ova_bright_med
ova_faint_84 = 3*(ova_faint_84-ova_faint_med) + ova_faint_med
ova_16 = 3*(ova_16-ova_med)+ ova_med
ova_bright_16 = 3*(ova_bright_16-ova_bright_med)+ ova_bright_med
ova_faint_16 = 3*(ova_faint_16-ova_faint_med)+ ova_faint_med
# plot as a function of radius
a = np.linspace(7.5,172.5,12)
#plt.rcParams['figure.figsize'] = [6, 4]
plt.plot(a,ova_med, color='tab:blue',label="Overdensity")
plt.plot(a,ova_bright_med,color= 'tab:red',label=r"Overdensity($\geq$5mJy)" )
plt.plot(a,ova_faint_med,color= 'tab:green',label="Overdensity(<5mJy)" )
plt.scatter(a,ova_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(a,ova_bright_med,color= 'tab:red')#,label="Overdensity(>4')" )
plt.scatter(a,ova_faint_med,color= 'tab:green')
plt.fill_between(a,ova_16,ova_84,color='tab:blue',alpha=0.1 )
plt.fill_between(a,ova_bright_16,ova_bright_84,color='tab:red',alpha=0.1 )
plt.fill_between(a,ova_faint_16,ova_faint_84,color='tab:green',alpha=0.1 )
#plt.tick_params(labelsize = 15)
plt.ylim(-0.9,4.4)
plt.legend()
plt.xlabel('Position angle (degree)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_a180.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_a180.eps',bbox_inches='tight') # traditional visualization
plt.close()


#plt.rcParams['figure.figsize'] = [6, 4]


plt.bar(a,ova_med,12,color='tab:blue',alpha = 0.6, label="Overdensity") # bar demonstration
plt.bar(a+3,ova_bright_med,6,color= 'tab:red',alpha = 0.6, label=r"Overdensity($\geq$5mJy)" )
plt.bar(a-3,ova_faint_med,6,color= 'tab:green',alpha = 0.6, label=r"Overdensity(<5mJy)" )
plt.errorbar(a,ova_med,yerr = [ova_med-ova_16,ova_84-ova_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '')
plt.errorbar(a+3,ova_bright_med,yerr =[ova_bright_med-ova_bright_16,ova_bright_84-ova_bright_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '' )
plt.errorbar(a-3,ova_faint_med,yerr =[ova_faint_med-ova_faint_16,ova_faint_84-ova_faint_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '' )
plt.ylim(-0.9,4.4)
#plt.tick_params(labelsize = 15)
plt.legend()
plt.xlabel('Position angle (degree)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_a_bar180.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_a_bar180.eps',bbox_inches='tight')
plt.close()   

