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
nca_out = np.zeros(12)
nca_in = np.zeros(12)
for i in range(12):
    for j in range(len(c_my)):
        if 15*i<=rg.position_angle(c_my[j]).degree<15*(i+1) or 180+15*i<=rg.position_angle(c_my[j]).degree<180+15*(i+1):
            nca[i] = nca[i]+1
            if c_my[j].separation(rg).arcmin<=4: 
                nca_in[i] = nca_in[i]+1
            else:
                nca_out[i] = nca_out[i]+1

na_sim = []
nain_sim = []
naout_sim = []
for i in range(ntrials):
    nca_sim = np.zeros(12)
    ncain_sim = np.zeros(12) # in sources > 5 mJy
    ncaout_sim = np.zeros(12)
    f_sim, c_sim = read_sim('mock_850/mock_map'+str(i)+'_rec.dat')
    for j in range(len(c_sim)):
        for k in range(12):
            if 15*k<=rg.position_angle(c_sim[j]).degree<15*(k+1) or 180+15*k<=rg.position_angle(c_sim[j]).degree<180+15*(k+1):
                nca_sim[k] = nca_sim[k]+1
                if c_sim[j].separation(rg).arcmin<=4:
                    ncain_sim[k] = ncain_sim[k]+1
                else:
                    ncaout_sim[k] = ncaout_sim[k]+1
    na_sim.append(nca_sim)
    nain_sim.append(ncain_sim)
    naout_sim.append(ncaout_sim)

na_sim = np.array(na_sim)
nain_sim = np.array(nain_sim)    
naout_sim = np.array(naout_sim) 
#np.save('angular/all.npy', np.array(na_sim))
np.save('angular/in.npy', np.array(nain_sim))
np.save('angular/out.npy', np.array(naout_sim))
# calculate the mean for every 20 maps
ova = []
ova_in = []
ova_out = []
ova_med = []
ova_in_med = []
ova_out_med = []
ova_84 = []
ova_in_84 = []
ova_out_84 = []
ova_16 = []
ova_in_16 = []
ova_out_16 = []
for i in range(int(ntrials/9)):
    na_sim_mean = np.nanmean(na_sim[i*9:(i+1)*9],axis=0) # each annulus mean value
    nain_sim_mean = np.nanmean(nain_sim[i*9:(i+1)*9],axis=0)
    naout_sim_mean = np.nanmean(naout_sim[i*9:(i+1)*9],axis=0)
    ova.append((nca-na_sim_mean)/na_sim_mean)
    ova_in.append((nca_in-nain_sim_mean)/nain_sim_mean)
    ova_out.append((nca_out-naout_sim_mean)/naout_sim_mean)

ova = np.array(ova)
ova_in = np.array(ova_in)
ova_out = np.array(ova_out)

for i in range(12):
    ova_med.append(np.nanmedian(ova[:,i]))
    ova_in_med.append(np.nanmedian(ova_in[:,i]))
    ova_out_med.append(np.nanmedian(ova_out[:,i]))
    ova_84.append(np.nanpercentile(ova[:,i],84.1))
    ova_in_84.append(np.nanpercentile(ova_in[:,i],84.1))
    ova_out_84.append(np.nanpercentile(ova_out[:,i],84.1))
    ova_16.append(np.nanpercentile(ova[:,i],15.9))
    ova_in_16.append(np.nanpercentile(ova_in[:,i],15.9))
    ova_out_16.append(np.nanpercentile(ova_out[:,i],15.9))

ova_med = np.array(ova_med)
ova_in_med = np.array(ova_in_med)
ova_out_med = np.array(ova_out_med)
ova_84 = np.array(ova_84)
ova_in_84 = np.array(ova_in_84)
ova_out_84 = np.array(ova_out_84)
ova_16 = np.array(ova_16)
ova_in_16 = np.array(ova_in_16)
ova_out_16 = np.array(ova_out_16)
ova_84 = 3*(ova_84-ova_med) + ova_med
ova_in_84 = 3*(ova_in_84-ova_in_med)+ ova_in_med
ova_out_84 = 3*(ova_out_84-ova_out_med) + ova_out_med
ova_16 = 3*(ova_16-ova_med)+ ova_med
ova_in_16 = 3*(ova_in_16-ova_in_med)+ ova_in_med
ova_out_16 = 3*(ova_out_16-ova_out_med)+ ova_out_med

# plot as a function of radius
a = np.linspace(7.5,172.5,12)
#plt.rcParams['figure.figsize'] = [6, 4]
plt.plot(a,ova_med, color='tab:blue',label="Overdensity")
plt.plot(a,ova_in_med,color= 'tab:green',label=r"Overdensity($\leq$4')" )
plt.plot(a,ova_out_med,color= 'tab:red',label="Overdensity(>4')" )
plt.scatter(a,ova_med, color='tab:blue')#,label="Overdensity(<4')")
plt.scatter(a,ova_in_med,color= 'tab:green')#,label="Overdensity(>4')" )
plt.scatter(a,ova_out_med,color= 'tab:red')
plt.fill_between(a,ova_16,ova_84,color='tab:blue',alpha=0.1 )
plt.fill_between(a,ova_in_16,ova_in_84,color='tab:green',alpha=0.1 )
plt.fill_between(a,ova_out_16,ova_out_84,color='tab:red',alpha=0.1 )
#plt.tick_params(labelsize = 15)
plt.legend()
plt.xlabel('Position angle (degree)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_a180r.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_a180r.eps',bbox_inches='tight') # traditional visualization
plt.close()


#plt.rcParams['figure.figsize'] = [6, 4]


plt.bar(a,ova_med,12,color='tab:blue',alpha = 0.6, label="Overdensity") # bar demonstration
plt.bar(a+3,ova_in_med,6,color= 'tab:green',alpha = 0.6, label=r"Overdensity($\leq$4')" )
plt.bar(a-3,ova_out_med,6,color= 'tab:red',alpha = 0.6, label="Overdensity(>4')"  )
plt.errorbar(a,ova_med,yerr = [ova_med-ova_16,ova_84-ova_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '')
plt.errorbar(a+3,ova_in_med,yerr =[ova_in_med-ova_in_16,ova_in_84-ova_in_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '' )
plt.errorbar(a-3,ova_out_med,yerr =[ova_out_med-ova_out_16,ova_out_84-ova_out_med],capsize=2,color = 'r',marker = 's',elinewidth = 1,markersize=3.5,linestyle = '' )
plt.ylim(-0.9,2.4)
#plt.tick_params(labelsize = 15)
plt.legend()
plt.xlabel('Position angle (degree)')#,fontsize=15)
plt.ylabel('Overdensity')#,fontsize=15)
plt.savefig('plot/ovd_a_bar180r.pdf',bbox_inches='tight')
plt.savefig('plot/ovd_a_bar180r.eps',bbox_inches='tight')
plt.close()   

