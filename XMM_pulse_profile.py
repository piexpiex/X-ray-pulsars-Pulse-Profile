import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.modeling import fitting,models
from stingray.events import EventList
from stingray.lightcurve import Lightcurve
from stingray.pulse.search import epoch_folding_search, z_n_search
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile
from pulse_profile import *

name_file="src_sd.events"
source='LMC X-4'
#datos (en principio)
asini=26.33 #[It-sec]
Porb=1.4084 #[days]
ecc=0.0
omega_d=0.0 #[degrees]
T0=51110.866 #[MJD]

period = 13.5

bin_time = 0.01

nbin = 40
nsinusoids=5

phase_profile=0.074195 #the phase obtained with NuSTAR_data_selector.py

####################
### Data lecture ###
####################

hdulist = fits.open(name_file)
hdulist.info()
header = hdulist[0].header

Tstart=hdulist['STDGTI04'].data['START']
Tstop=hdulist['STDGTI04'].data['STOP']

#PI
hdulist = fits.open("XMMpha.fits")
PI = hdulist[0].data
#times
hdulist = fits.open("XMM.fits")
times = hdulist[0].data
T_star_stop=[[Tstart[0],Tstop[0]]]
Total_exptime=0
for j in range(len(Tstart)):
	Total_exptime=Total_exptime+Tstop[j]-Tstart[j]
	if j>0:
		T_star_stop.append([Tstart[j],Tstop[j]])
T_star_stop=np.array(T_star_stop)

#####################
### Pulse profile ###
#####################

PI_ranges=[95,600,1200,2400]
Energy_ranges=[0.5,3,6,12] #KeV

#PHA 95-600
A1=np.where((PI<600) & (PI>95))
E1=times[A1]
#PHA 600-1200
A2=np.where((PI<1200) & (PI>600))
E2=times[A2]
#PHA 1200-2400
A3=np.where((PI<2400) & (PI>1200))
E3=times[A3]

#para las unidades reales
#Err=profile_err*nbin/Total_exptime
#profile=profile*nbin/Total_exptime

Pulse_total=pulse_profile()
Pulse_total.times=times
Pulse_total.phase_profile=phase_profile
Pulse_total.nbin=nbin
Pulse_total.T_star_stop=T_star_stop
Pulse_total.ph, Pulse_total.profile, Pulse_total.profile_err = Pulse_total.profile()

_ = plot_profile(Pulse_total.ph, Pulse_total.profile,Pulse_total.profile_err)
plt.show()

Pulse_total.ph, Pulse_total.profilenorm, Pulse_total.profile_err = Pulse_total.profile_norm()
ph2, ffit,A,F,Sigma = Pulse_total.adjusment(5)
print('ppppl',A,F,Sigma)

#1
Pulse_total_1=pulse_profile()
Pulse_total_1.times=E1
Pulse_total_1.phase_profile=phase_profile
Pulse_total_1.nbin=nbin
Pulse_total_1.T_star_stop=T_star_stop
Pulse_total_1.ph, Pulse_total_1.profile, Pulse_total_1.profile_err = Pulse_total_1.profile()
Pulse_total_1.ph, Pulse_total_1.profilenorm, Pulse_total_1.profile_err = Pulse_total_1.profile_norm()
ph2, ffit_1,A_1,F_1,Sigma_1 = Pulse_total_1.adjusment(5)
#2
Pulse_total_2=pulse_profile()
Pulse_total_2.times=E2
Pulse_total_2.phase_profile=phase_profile
Pulse_total_2.nbin=nbin
Pulse_total_2.T_star_stop=T_star_stop
Pulse_total_2.ph, Pulse_total_2.profile, Pulse_total_2.profile_err = Pulse_total_2.profile()
Pulse_total_2.ph, Pulse_total_2.profilenorm, Pulse_total_2.profile_err = Pulse_total_2.profile_norm()
ph2, ffit_2,A_2,F_2,Sigma_2 = Pulse_total_2.adjusment(5)
#3
Pulse_total_3=pulse_profile()
Pulse_total_3.times=E3
Pulse_total_3.phase_profile=phase_profile
Pulse_total_3.nbin=nbin
Pulse_total_3.T_star_stop=T_star_stop
Pulse_total_3.ph, Pulse_total_3.profile, Pulse_total_3.profile_err = Pulse_total_3.profile()
Pulse_total_3.ph, Pulse_total_3.profilenorm, Pulse_total_3.profile_err = Pulse_total_3.profile_norm()
ph2, ffit_3,A_3,F_3,Sigma_3 = Pulse_total_3.adjusment(5)

plt.subplot(2,2,1)



plt.suptitle('Source:'+source+'  \n XMM-Newton observations \n Pulse period '+str(period)+'s',fontsize=12)
plt.step(Pulse_total.ph,Pulse_total.profilenorm,where='mid',color='k',label='total pulse')
plt.errorbar(Pulse_total.ph,Pulse_total.profilenorm,yerr=Pulse_total.profile_err,fmt='ko',markersize=0.5)
plt.ylim(-15,15)
plt.plot(ph2,ffit(ph2), 'k--')
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)


plt.subplot(2,2,2)
plt.step(Pulse_total_1.ph,Pulse_total_1.profilenorm,where='mid',color='y',label='0.5-3 KeV')
plt.errorbar(Pulse_total_1.ph,Pulse_total_1.profilenorm,yerr=Pulse_total_1.profile_err,fmt='yo',markersize=0.5)

plt.plot(ph2,ffit_1(ph2), 'y--')
plt.ylim(-15,15)
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,2,3)
plt.step(Pulse_total_2.ph,Pulse_total_2.profilenorm,where='mid',color='b',label='3-6 KeV')
plt.errorbar(Pulse_total_2.ph,Pulse_total_2.profilenorm,yerr=Pulse_total_2.profile_err,fmt='bo',markersize=0.5)

plt.plot(ph2,ffit_2(ph2), 'b--')
plt.ylim(-15,15)
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,2,4)
plt.step(Pulse_total_3.ph,Pulse_total_3.profilenorm,where='mid',color='g',label='6-12 KeV')
plt.errorbar(Pulse_total_3.ph,Pulse_total_3.profilenorm,yerr=Pulse_total_3.profile_err,fmt='go',markersize=0.5)

plt.ylim(-15,15)
plt.plot(ph2,ffit_3(ph2), 'g--')
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
plt.show()
#saving the parameters

#A_3,F_3,Sigma_3
AD1=np.concatenate((A_1,F_1,Sigma_1),axis=None)
AD2=np.concatenate((A_2,F_2,Sigma_2),axis=None)
AD3=np.concatenate((A_3,F_3,Sigma_3),axis=None)
print('popop')
print(AD1)
print(AD2)
print(AD3)
Parameters_XMM=np.concatenate((AD1,AD2,AD3),axis=None)

Parameters_nustar=np.array(Parameters_XMM)
print(Parameters_XMM)
hdu = fits.PrimaryHDU(Parameters_XMM)
hdul = fits.HDUList([hdu])
hdul.writeto('Parameters_XMM.fits',overwrite=True)
