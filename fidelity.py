import numpy as np

path = input('Path of mocked maps (e.g. mock_850/): ')
nu = int(input('Number of mocked maps: '))
#step = int(input('The step size of the FWHM: '))
#num = int(input('Number of different kernels: '))
#ini = int(input('Initial kernel size: '))
flux_min = float(input('Flux min: '))
flux_max = float(input('Flux max: '))
noise_min = float(input('Noise min: '))
noise_max = float(input('Noise max: '))

flux = np.linspace(flux_min,flux_max,26) 
snr = np.arange(0,26,1)
noise = np.linspace(noise_min,noise_max,26) 

#n_recflux = []
#n_spuflux = []
#n_recsnr = []
#n_spusnr = []
#n_recnoise = []
#n_spunoise = []

#for j in range(num):
flux_rec = np.zeros(25)
flux_spu = np.zeros(25)
snr_rec = np.zeros(25)
snr_spu = np.zeros(25)
noise_rec = np.zeros(25)
noise_spu = np.zeros(25)
for i in range(nu):
    try:
        file_rec = path+ 'mock_map'+str(i)+'_rec.dat'
        file_spu = path+'spu_map'+str(i)+'.dat'
        rec = np.loadtxt(file_rec,skiprows=1)
        spu = np.loadtxt(file_spu,skiprows=1)     
        if rec.ndim >1 and rec.shape[0] !=0:
            rec_flux = rec[:,2]
            rec_noise = rec[:,3]
            rec_snr = rec[:,4]
            for k in range(25):
                for l in range(rec.shape[0]):
                    if flux[k]<=rec_flux[l]<flux[k+1]:
                        flux_rec[k]=flux_rec[k]+1
                    if snr[k]<=rec_snr[l]<snr[k+1]:
                        snr_rec[k]=snr_rec[k]+1
                    if noise[k]<=rec_noise[l]<noise[k+1]:
                        noise_rec[k]=noise_rec[k]+1
        elif rec.ndim ==1 and rec.shape[0] !=0:
            rec_flux = rec[2]
            rec_noise = rec[3]
            rec_snr = rec[4]
            for k in range(25):
                if flux[k]<=rec_flux<flux[k+1]:
                    flux_rec[k]=flux_rec[k]+1
                if snr[k]<=rec_snr<snr[k+1]:
                    snr_rec[k]=snr_rec[k]+1
                if noise[k]<=rec_noise <noise[k+1]:
                    noise_rec[k]=noise_rec[k]+1
        if spu.ndim >1 and spu.shape[0] !=0:
            spu_flux = spu[:,2]
            spu_noise = spu[:,3]
            spu_snr = spu[:,4]
            for k in range(25):
                for l in range(spu.shape[0]):
                    if flux[k]<=spu_flux[l]<flux[k+1]:
                        flux_spu[k]=flux_spu[k]+1
                    if snr[k]<=spu_snr[l]<snr[k+1]:
                        snr_spu[k]=snr_spu[k]+1
                    if noise[k]<=spu_noise[l]<noise[k+1]:
                        noise_spu[k]=noise_spu[k]+1
        elif spu.ndim ==1 and spu.shape[0] !=0:
            spu_flux = spu[2]
            spu_noise = spu[3]
            spu_snr = spu[4]
            for k in range(25):
                if flux[k]<=spu_flux<flux[k+1]:
                    flux_spu[k]=flux_spu[k]+1
                if snr[k]<=spu_snr<snr[k+1]:
                    snr_spu[k]=snr_spu[k]+1
                if noise[k]<=spu_noise <noise[k+1]:
                    noise_spu[k]=noise_spu[k]+1
    except:
        print(str(i)+' map not exist!')
n_recflux = flux_rec#.append(flux_rec)
n_spuflux = flux_spu#.append(flux_spu)
n_recsnr = snr_rec#.append(snr_rec)
n_spusnr = snr_spu#.append(snr_spu)
n_recnoise = noise_rec#.append(noise_rec)
n_spunoise = noise_spu#.append(noise_spu)    
'''
ff = []
fs = []
fn = []
for j in range(num):
'''
fide_f = np.zeros(25)
fide_s = np.zeros(25)
fide_n = np.zeros(25)

for i in range(25):
    if n_recflux[i]!=0:
        fide_f[i] = 1-n_spuflux[i]/n_recflux[i]
    else:
        fide_f[i] = np.nan
    if n_recsnr[i]!=0:
        fide_s[i] = 1-n_spusnr[i]/n_recsnr[i]
    else:
        fide_s[i] = np.nan
    if n_recnoise[i]!=0:
        fide_n[i] = 1-n_spunoise[i]/n_recnoise[i]
    else:
        fide_n[i] = np.nan
ff = (flux[:-1],fide_f) #.append((flux[:-1],fide_f))
fs = (snr[:-1],fide_s) #.append((snr[:-1],fide_s))
fn = (noise[:-1],fide_n) #.append((noise[:-1],fide_n))
        
np.save('fide_snr'+str(path[4:-1])+'.npy', np.array(fs))
np.save('fide_flux'+str(path[4:-1])+'.npy', np.array(ff))
np.save('fide_noise'+str(path[4:-1])+'.npy', np.array(fn))

