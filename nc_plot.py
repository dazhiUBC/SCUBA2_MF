from sim_proto import *
import pandas as pd 


'''not sure if it is finished, probably go to the notebook'''
def nc_reader(num):
    file = 'sim/Ture_NC_'+str(num)+'deboosted_sources_4C23_850_cal_crop_MF.dat' 
    nc = pd.read_csv(file,delimiter=' ')
    return nc['cor_num'],nc['cor_num_in']
'''
def nc_reader9(num):
    file = '/sim/Ture_NC_'+str(num)+'deboosted_sources_4C23_850_cal_crop_MF.dat' 
    nc = pd.read_csv(file,delimiter=' ')
    return nc['cor_num']
'''
def ncs(s, n0=7180, s0=2.5, gamma=1.5):
    '''scuba NC by Geach+17'''
    return (n0/s0)*(s/s0)**(-gamma) *np.exp(-s/s0)

def nca(s,a=5.9,b=0.4,s0=5.1,n0=1200):
    '''ALMA NC by Stach+18'''
    return (n0/s0)*((s/s0)**a+(s/s0)**b)**(-1)


def nc9(num):
    '''normalize the number count whole region'''
    area = np.pi * (9/60)**2
    return num/area

def nc4(num):
    '''normalize the number count inner region'''
    area = np.pi * (4/60)**2
    return num/area

inn = []
out = []
for i in range(10000):
    try:
        out_i, inn_i = nc_reader(i)
        inn.append(inn_i)
        out.append(out_i)
    except:
        True
NC = np.array(out)
NC_in = np.array(inn)


nc_mean = [] # overall number count 
std = []
nc4mean = []
std4 = []
for i in range(15):
    ncc.append(np.mean(NC[:,i]))
    std.append(np.std(NC[:,i]))
    nc4mean.append(np.mean(NC_in[:,i+3]))
    std4.append(np.std(NC_in[:,i+3]))


    
me = []
me_std = []
for i in range(15):
    me.append(nc(nc4mean[i]))
    me_std.append(nc(std4[i]))
    
    
me9 = []
me_std9 = []
for i in range(15):
    me9.append(nc9(ncc[i]))
    me_std9.append(nc9(std[i]))
    
me = np.array(me)
me_std = np.array(me_std)

me9 = np.array(me9)
me_std9 = np.array(me_std9)

s = np.linspace(1,20,20)

n = nc(s)
nbc = ncd(s)
e = np.linspace(1.5,15.5,15)
plt.errorbar(e[:-1],me[:-1],yerr=0.5*me_std[:-1],fmt='o',color='b',ecolor='b',mew=1,capsize=4,label=r"Number counts for 4C23.56 (r<3')")
plt.errorbar(e[:-1],me9[:-1],yerr=0.8*me_std9[:-1],fmt='o',color='r',ecolor='r',mew=1,capsize=4,label=r"Number counts for 4C23.56 (r<9')")
plt.plot(s,n,'tab:orange',label=r'SCUBA-2 850$\mu m$ (Geach+2017)')
plt.plot(s,nbc,'tab:blue',label=r'ALMA 870$\mu m$ (Stach+2018)')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1.5,50)
plt.ylabel(r'dN/dS (degree$^{-2}$ mJy$^{-1}$)')
plt.xlabel(r'S$_{850}$ (mJy)')
plt.ylim(2,3000)
plt.legend()
plt.savefig('nc.pdf',bbox_inches='tight')
plt.savefig('nc.eps',bbox_inches='tight')