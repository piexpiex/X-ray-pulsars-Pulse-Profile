import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from stingray.events import EventList
#from stingray.lightcurve import Lightcurve
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile
from binary_cor import *
from read_files import *
import os
import warnings
warnings.filterwarnings('ignore')

READ=read_files()

name_file=READ[1] #name of the files with the data

#datos (en principio)
asini=READ[3] #[It-sec]
Porb=READ[4] #[days]
ecc=READ[5]
omega_d=READ[6] #[degrees]
T0=READ[7] #[MJD]

period = READ[8]

nbin = READ[12]

overwrite=READ[14]
Z_2_check=READ[15]

times=np.array([])

key_overwrite=0

if overwrite=='n' and READ[10]!=0:
	try:
		fits.open('fits_folder/XMM_pha.fits')
		fits.open('fits_folder/XMM_times.fits')
		key_overwrite=1
	except:
		print('XMM files not found')
		print('I proceed to find the exact frequencies')
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

PHA=np.array([])
print(name_file)
for k in range(len(name_file)):
	print('file',k+1,'of',len(name_file))
	hdulist = fits.open(name_file[k])
	hdulist.info()
	header = hdulist[0].header
	DATA = hdulist[1].data
	times_k=np.zeros((len(DATA)))
	PHA_k=np.zeros((len(DATA)))
	Column=np.where(np.array(DATA.columns.names)=='PHA')
	TIME_column=Column[0][0]
	Column_2=np.where(np.array(DATA.columns.names)=='PHA')
	PHA_column=Column_2[0][0]
	for j in range(len(times_k)):
		if j==round(len(DATA)*0.05):print(5,"%")
		if j==round(len(DATA)*0.25):print(25,"%")
		if j==round(len(DATA)*0.5):print(50,"%")
		if j==round(len(DATA)*0.8):print(80,"%")
		if j==len(DATA)-1:print(100,"%")
		DATAj=DATA[j]
		times_k[j]=DATAj[0]#[TIME_column]
		PHA_k[j]=DATAj[PHA_column]
	times=np.append(times,times_k)
	PHA=np.append(PHA,PHA_k)

hdu = fits.PrimaryHDU(PHA)
hdul = fits.HDUList([hdu])
hdul.writeto('fits_folder/XMM_pha.fits',overwrite=True)
print('Correction of the binary sistem delay')
TIME=Binary_orbit(time=times,asini=asini,ecc=ecc,porb=Porb,omega_d=omega_d ,t0=T0)
hdu = fits.PrimaryHDU(TIME)
hdul = fits.HDUList([hdu])
hdul.writeto('fits_folder/XMM_times.fits',overwrite=True)

times=TIME

# We will search for pulsations over a range of frequencies around the known pulsation period.
obs_length = times[len(times)-1]-times[0]
df_min = 10/obs_length
oversampling=15
df = df_min / oversampling
frequencies = np.arange(1/period - 500 * df, 1/period + 500 * df, df)

freq, efstat = epoch_folding_search(times, frequencies, nbin=nbin)

pulse_frequency=freq[np.where(efstat==max(efstat))[0][0]]
duplication=1
while max(efstat)<200:
	pulse_frequency_value=input('the value of the pulse frequency is small, do you want to search a better value (Y/N)')
	if pulse_frequency_value=='Y' or pulse_frequency_value=='y':
		duplication=duplication*2
		frequencies = np.append(np.arange(1/period - 200*duplication * df,1/period - 200*duplication/2 * df,df), np.arange(1/period + 200*duplication/2 * df,1/period + 200*duplication * df, df))
		freq, efstat = epoch_folding_search(times, frequencies, nbin=nbin)
		pulse_frequency=freq[np.where(efstat==max(efstat))[0][0]]
	elif pulse_frequency_value=='N' or pulse_frequency_value=='n':
		break
	else:
		pulse_frequency_value=input('please Y or N')
print('pulse frequency',pulse_frequency)
_=write_files('pulse_frequency_XMM',pulse_frequency)
# ---- PLOTTING --------
plt.figure()
plt.plot(freq, efstat, label='EF statistics')
plt.axhline(nbin - 1, ls='--', lw=3, color='k', label='n - 1')
plt.xlabel('Frequency (Hz)')
plt.ylabel('EF Statistics')
_ = plt.legend()
plt.savefig("figures_folder/XMM_pulse_frequency_search.pdf")

if Z_2_check=='N':
	exit()

# We will search for pulsations over a range of frequencies around the known pulsation period.
nharm = 1
freq, zstat = z_n_search(times, frequencies, nbin=nbin, nharm=nharm)

# ---- PLOTTING --------
plt.figure()
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
plt.savefig("figures_folder/XMM_pulse_frequency_search_check.pdf")
