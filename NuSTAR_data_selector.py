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
from fpropias import *
from scipy.optimize import curve_fit

name_file=["nu30102041002A01_cl_src_bary.evt","nu30102041002B01_cl_src_bary.evt"] #name of the file with the data

hdulist = fits.open(name_file[0])
hdulist.info()
header = hdulist[2].header
print('...')
#print(header)
print('...')

#datos (en principio)
asini=26.33 #[It-sec]
Porb=1.4084 #[days]
ecc=0.0
omega_d=0.0 #[degrees]
T0=51110.866 #[MJD]

mpl.rcParams['figure.figsize'] = (10, 6)
def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))

period = 13.5

bin_time = 0.01

nbin = 40

hdulist = fits.open(name_file[0])
hdulist.info()
header = hdulist[0].header

#lo anterior
ev =EventList()
ev = ev.read(name_file[0], 'fits')
times=ev.time

#######################################
if len(name_file)>1:
	hdulist = fits.open(name_file[1])
	hdulist.info()
	header = hdulist[2].header

	hdulist = fits.open(name_file[1])
	hdulist.info()
	header = hdulist[0].header

	#lo anterior
	ev =EventList()	
	ev = ev.read(name_file[1], 'fits')

	times2=ev.time
#######################################
times=np.append(times,times2)

PI=ev.pi
CT=ev.ncounts
Tstart=hdulist['GTI'].data['START']
Tstop=hdulist['GTI'].data['STOP']
TIME=Binary_orbit(time=times,asini=asini,ecc=ecc,porb=Porb,omega_d=omega_d ,t0=T0)
hdu = fits.PrimaryHDU(TIME)
hdul = fits.HDUList([hdu])
hdul.writeto('nustar.fits')

# We will search for pulsations over a range of frequencies around the known pulsation period.
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
