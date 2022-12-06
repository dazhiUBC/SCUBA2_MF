from completeness import * 
import matplotlib.pyplot as plt

path = input('The name of the directory (e.g. mock_850): ') + '/'
n_trials = int(input('The number of mock maps: '))

def match_dist(inp_file, rec_file, completeness=False):
    '''match the input source with the recovered sources and return the matched catalog with distance'''
    inp_ra, inp_dec, inp_flux, inp_noise, rec_ra, rec_dec, rec_flux, rec_err = get_inputs(inp_file, rec_file)
    inp_ra = list(inp_ra)
    inp_dec = list(inp_dec)
    inp_flux = list(inp_flux)
    inp_noise = list(inp_noise)

    data = zip(rec_ra, rec_dec, rec_flux)
    data2 = sorted(data, key=lambda tup: tup[2], reverse=True)
    rec_ra, rec_dec, rec_flux = zip(*data2)
    rec_coord = SkyCoord(np.array(rec_ra)*u.degree, np.array(rec_dec)*u.degree)

    '''loops through outputs'''
    for r_c, r_f in zip(rec_coord, rec_flux):
        '''searches inputs to find match to recovered'''
        inp_coord = SkyCoord(np.array(inp_ra)*u.degree, np.array(inp_dec)*u.degree)
        d2d = r_c.separation(inp_coord).arcsecond

        matches = []
        ind = 0
        for d in d2d:
          if d < 5:
            i_f = inp_flux[ind]
            i_noise = inp_noise[ind]
            matches.append((i_f, i_noise,d, ind))
          ind+=1

        if len(matches) > 0:
          match = max(matches, key = lambda t: t[0])
          catalog.append((match[0], match[1], r_f, match[2]))

          '''remove the matched value (avoid doubles)'''
          index = match[3]
          del inp_ra[index]
          del inp_dec[index]
          del inp_flux[index]
          del inp_noise[index]
    return catalog
      


catalog = []


for i in range(n_trials): 
    inp_file = path+'mock_map'+str(i)+'.dat'
    rec_file =path+'mock_map'+str(i)+'_rec.dat'
    match_dist(inp_file, rec_file, completeness=True)
    
snr = np.array(catalog)[:,2]/np.array(catalog)[:,1]
offset = np.array(catalog)[:,3]

off_bin = []
for j in range(10):
    snr_bin = []
    for i in range(len(catalog)):
        if 3.5+j<=snr[i]<4.5+j:
            snr_bin.append(offset[i])
    off_bin.append(snr_bin)
    
mean = []
med = []
o84 = []
o16 = []
for i in range(10):
    mean.append(np.mean(off_bin[i]))
    med.append(np.median(off_bin[i]))
    o84.append(np.percentile(off_bin[i],84.1))            
    o16.append(np.percentile(off_bin[i],15.9))
    
x = np.linspace(4,13,10)
plt.plot(x, med, color = 'r')
plt.fill_between(x,o16,o84,color='r',alpha=0.1,label=r'1$\sigma$ deviation')
plt.scatter(snr,offset,s=0.05,alpha=0.3)

plt.ylabel('Offset (arcsecond)')
plt.xlabel('SNR')
plt.ylim(-0.4,5)
plt.xlim(4,13)
plt.legend()
plt.savefig('offset.pdf',bbox_inches='tight')