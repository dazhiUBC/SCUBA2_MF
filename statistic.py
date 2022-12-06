import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fid_snr = np.load('fide_snr_850.npy')
de = pd.read_csv('850/deboosting_values_snr.dat',delimiter = ' ')
com = pd.read_csv('850/completeness_values_snr.dat',delimiter = ' ')
plt.plot(com['input_snr'],com['fraction'],color='tab:green',label = 'Completeness')
plt.plot(de['observed_snr'], de['deboost_value_mean'],color='tab:red',label = 'Deboosting')
plt.fill_between(de['observed_snr'],de['deboost_value_mean']+de['std'],de['deboost_value_mean']-de['std'],color='tab:red',alpha=0.1)
plt.plot(fid_snr[0],fid_snr[1],label = 'Fidelity',color='tab:blue')
plt.xlim(3.5,16)
plt.ylim(0.31,1.08)
plt.xlabel(r'SNR$_{rec}$')
plt.ylabel('Factor')
plt.legend()
plt.savefig('deboost_new.pdf',bbox_inches='tight')