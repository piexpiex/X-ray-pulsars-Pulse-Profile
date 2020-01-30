import numpy as np
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from matplotlib import pyplot as plt
from astropy.io import fits
import seaborn as sb
import matplotlib as mpl
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile
from binary_cor import *
from scipy.optimize import curve_fit
from math import isnan

name_file="src_sd.events"



#datos (en principio)
asini=26.33 #[It-sec]
Porb=1.4084 #[days]
ecc=0.0
omega_d=0.0 #[degrees]
T0=51110.866 #(MJD)

period = 13.5

bin_time = 0.01

nbin = 20

hdulist = fits.open(name_file)
hdulist.info()
header = hdulist[0].header


DATA = hdulist[1].data
print(DATA.columns.names)
times=np.zeros((len(DATA)))
PHA=np.zeros((len(DATA)))
nmn=np.where(np.array(DATA.columns.names)=='PHA')
PHA_column=nmn[0][0]
print(PHA_column)

print('------------------------------------')
for j in range(len(times)):
	if j==round(len(DATA)*0.05):print(5,"%")
	if j==round(len(DATA)*0.25):print(25,"%")
	if j==round(len(DATA)*0.5):print(50,"%")
	if j==round(len(DATA)*0.8):print(80,"%")
	if j==len(DATA)-1:print(100,"%")
	DATAj=DATA[j]
	times[j]=DATAj[0]
	PHA[j]=DATAj[PHA_column]

hdu = fits.PrimaryHDU(PHA)
hdul = fits.HDUList([hdu])
hdul.writeto('XMMpha.fits',overwrite=True)

for j in range(len(times)):
	if isnan(times[j])==True:
		print('DANGER --> lecture error at line',j,'of the file')

TIME=Binary_orbit(time=times,asini=asini,ecc=ecc,porb=Porb,omega_d=omega_d ,t0=T0)
hdu = fits.PrimaryHDU(TIME)
hdul = fits.HDUList([hdu])
hdul.writeto('XMM.fits',overwrite=True)

times=TIME

# We will search for pulsations over a range of frequencies around the known pulsation period.
obs_length = times[len(times)-1]-times[0]
df_min = 1/obs_length
oversampling=15
df = df_min / oversampling
frequencies = np.arange(1/period - 200 * df, 1/period + 200 * df, df)

freq, efstat = epoch_folding_search(times, frequencies, nbin=nbin)

# ---- PLOTTING --------
plt.figure()
plt.plot(freq, efstat, label='EF statistics')
plt.axhline(nbin - 1, ls='--', lw=3, color='k', label='n - 1')
plt.axvline(1/period, lw=3, alpha=0.5, color='r', label='Correct frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('EF Statistics')
_ = plt.legend()
plt.show()


# We will search for pulsations over a range of frequencies around the known pulsation period.
nharm = 1
freq, zstat = z_n_search(times, frequencies, nbin=nbin, nharm=nharm)

# ---- PLOTTING --------
plt.figure()
plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

plt.axvline(1/period, color='r', lw=3, alpha=0.5, label='Correct frequency')
plt.xlim([frequencies[0], frequencies[-1]])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Statistics - d.o.f.')
plt.legend()
plt.figure(figsize=(15, 5))
plt.plot(freq, (zstat - nharm), label='$Z_2$ statistics')
plt.plot(freq, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

plt.axvline(1/period, color='r', lw=3, alpha=0.5, label='Correct frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Statistics - d.o.f. (Zoom)')

plt.ylim([-15, 15])
_ = plt.xlim([frequencies[0], frequencies[-1]])
plt.show()
