import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
from stingray.events import EventList
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile
from binary_cor import *
from read_files import *
import os
import warnings
warnings.filterwarnings('ignore')

READ=read_files() # run.sh data

name_file=READ[1] 
if name_file[0]==' ':
	exit()
source=READ[2]
asini=READ[3] 
Porb=READ[4]
ecc=READ[5]
omega_d=READ[6] 
T0=READ[7]

period = READ[8]

nbin = READ[12]

overwrite=READ[14]
Z_2_check=READ[15]

period_ranges=READ[16]
period_ranges=np.array(period_ranges)
period_ranges=period_ranges.astype(np.float)

period_bins=READ[17]
XMM_time=READ[18]

key_overwrite=0

if overwrite=='n' and READ[10]!=0:
	try:
		fits.open('fits_folder/'+add_space(source)+'_XMM_times.fits')
		key_overwrite=1
	except:
		print('XMM files not found')
		key_overwrite=0
		
if key_overwrite==1:
	exit()

try:
	os.mkdir('fits_folder')
	os.mkdir('figures_folder')
	print('making directories')
	print('starting the analysis')
except:
	print('starting the analysis')
	
times=np.array([])
PHA=np.array([])

for k in range(len(name_file)):
	print('file',k+1,'of',len(name_file))
	hdulist = fits.open(name_file[k])
	hdulist.info()
	header = hdulist[0].header
	data = hdulist[1].data
	Column=np.where(np.array(data.columns.names)=='TIME')
	TIME_column=Column[0][0]
	Column_2=np.where(np.array(data.columns.names)=='PHA')
	PHA_column=Column_2[0][0]
	times_k=data.field(TIME_column)
	PHA_k=data.field(PHA_column) 
	times=np.append(times,times_k)
	PHA=np.append(PHA,PHA_k)

if asini==0 and Porb==0 and ecc==0 and omega_d==0 and T0==0:
	print('No possible correction of the binary sistem delay')
	TIME=times
else:
	print('Correction of the binary sistem delay')
	TIME=Binary_orbit(time=times,asini=asini,ecc=ecc,porb=Porb,omega_d=omega_d ,t0=T0)

c1 = fits.Column(name='TIMES', array=TIME, format='D')
c2 = fits.Column(name='PHA',array=PHA, format='D')

t = fits.BinTableHDU.from_columns([c1,c2],name='VALUES')
t.writeto('fits_folder/'+add_space(source)+'_XMM_times.fits',overwrite=True)

times=TIME

print('Search for the best frequency')

df=(period_ranges[1]-period_ranges[0])/period_bins
frequencies = 1/np.arange(period_ranges[0], period_ranges[1], df)

freq, efstat = epoch_folding_search(times, frequencies, nbin=nbin)

pulse_frequency=freq[np.where(efstat==max(efstat))[0][0]]

print('pulse frequency',pulse_frequency)
_=write_files('pulse_frequency_XMM',pulse_frequency)

#fitting epoch folding distribution with Lorentzian curve
g_init=models.Lorentz1D(amplitude=max(efstat)-min(efstat),x_0=pulse_frequency,fwhm=pulse_frequency/500)+models.Const1D(amplitude=min(efstat))
fit_g=fitting.LevMarLSQFitter()
bin_max=[np.where(efstat==max(efstat))[0][0]][0]
bin_left_min=0
bin_right_min=len(efstat)


for j in range(min(bin_max,len(efstat)-bin_max)-2):
	d_efstat=efstat[bin_max-j-1]-efstat[bin_max-j]
	if d_efstat>0 and max(efstat)-efstat[bin_max-j]>(max(efstat)-min(efstat))/5:
		bin_left_min=bin_max-j
		break
for j in range(min(bin_max,len(efstat)-bin_max)-2):
	d_efstat=efstat[bin_max+j+1]-efstat[bin_max+j]
	if d_efstat>0 and max(efstat)-efstat[bin_max+j]>(max(efstat)-min(efstat))/5:
		bin_right_min=bin_max+j
		break	

# ---- PLOTTING --------
plt.figure()
plt.suptitle('Source:'+source+'  \n XMM-Newton observations',fontsize=12)
plt.plot(freq, efstat, label='EF statistics')
try:
	g=fit_g(g_init,freq[bin_left_min:bin_right_min+1],efstat[bin_left_min:bin_right_min+1])
	#print(g[0].amplitude[0],g[0].stddev[0],g[0].mean[0],g[1].amplitude[0])
	fwhm=g[0].fwhm[0]#*2.35
	print('pulse frequency error',fwhm/2)
	plt.plot(freq[bin_left_min:bin_right_min],g(freq[bin_left_min:bin_right_min]),'b',label='Pulse frequency='+str(round(pulse_frequency ,2+int(-np.log10(fwhm/2))))+' $\pm$ '+str(round(fwhm/2,2+int(-np.log10(fwhm/2)))))
except:
	pass
plt.axhline(nbin - 1, ls='--', lw=3, color='k', label='n - 1')
plt.xlabel('Frequency (Hz)')
plt.ylabel('EF Statistics')
_ = plt.legend()
plt.savefig('figures_folder/'+add_space(source)+'_XMM_pulse_frequency_search.pdf')

if Z_2_check=='N':
	exit()

# We will search for pulsations over a range of frequencies around the known pulsation period.
nharm = 1
freq, zstat = z_n_search(times, frequencies, nbin=nbin, nharm=nharm)

# ---- PLOTTING --------
plt.figure()
plt.suptitle('Source:'+source+'  \n XMM-Newton observations',fontsize=12)
plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

plt.xlim([frequencies[0], frequencies[-1]])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Statistics - d.o.f.')
plt.legend()
plt.figure(figsize=(15, 5))
plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

plt.xlabel('Frequency (Hz)')
plt.ylabel('Statistics - d.o.f. (Zoom)')

plt.ylim([-15, 15])
_ = plt.xlim([frequencies[0], frequencies[-1]])
plt.savefig('figures_folder/'+add_space(source)+'_XMM_pulse_frequency_search_check.pdf')
