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
from read_files import *

READ=read_files()
name_file=READ[1] #name of the files with the data

#datos (en principio)
source=READ[2]
asini=READ[3] #[It-sec]
Porb=READ[4] #[days]
ecc=READ[5]
omega_d=READ[6] #[degrees]
T0=READ[7] #[MJD]

period = READ[8]

bin_time = READ[9]

nbin = READ[10]

nsinusoids=READ[11]

phase_profile=READ[13] #the phase obtained with NuSTAR_data_selector.py

####################
### Data lecture ###
####################

PI=np.array([])

for j in range(len(name_file)):
	hdulist = fits.open(name_file[j])
	hdulist.info()
	header = hdulist[0].header

	#lo anterior
	ev =EventList()	
	ev = ev.read(name_file[j], 'fits')
	PI=np.append(PI,ev.pi)

Tstart=hdulist['GTI'].data['START']
Tstop=hdulist['GTI'].data['STOP']
T_star_stop=[[Tstart[0],Tstop[0]]]
Total_exptime=0
for j in range(len(Tstart)):
	Total_exptime=Total_exptime+Tstop[j]-Tstart[j]
	if j>0:
		T_star_stop.append([Tstart[j],Tstop[j]])
T_star_stop=np.array(T_star_stop)

#times
hdulist = fits.open("nustar.fits")
times = hdulist[0].data

#####################
### Pulse profile ###
#####################

#PHA 35-110
#A1=np.where((PI<110) & (PI>35))
#E1=times[A1]
#PHA 110-260
#A2=np.where((PI<260) & (PI>110))
#E2=times[A2]
#PHA 260-650
#A3=np.where((PI<650) & (PI>260))
#E3=times[A3]
#PHA 650-1910
#A4=np.where((PI<1910) & (PI>650))
#E4=times[A4]

Energy_ranges=[3,4.5,6,9,12,15,18,27.6,78] #KeV

def PI_calculator(E):
	PI_range=[0.0]*len(E)
	for j in range(len(E)):
			PI_range[j]= -40+25*E[j]
	return(PI_range)

PI_ranges=PI_calculator(Energy_ranges)

A1=np.where((PI<72.5) & (PI>35))
E1=times[A1]
A2=np.where((PI<110) & (PI>72.5))
E2=times[A2]
A3=np.where((PI<185) & (PI>110))
E3=times[A3]
A4=np.where((PI<260) & (PI>185))
E4=times[A4]
A5=np.where((PI<335) & (PI>260))
E5=times[A5]
A6=np.where((PI<410) & (PI>335))
E6=times[A6]
A7=np.where((PI<650) & (PI>410))
E7=times[A7]
A8=np.where((PI<1910) & (PI>650))
E8=times[A8]
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
#4
Pulse_total_4=pulse_profile()
Pulse_total_4.times=E4
Pulse_total_4.phase_profile=phase_profile
Pulse_total_4.nbin=nbin
Pulse_total_4.T_star_stop=T_star_stop
Pulse_total_4.ph, Pulse_total_4.profile, Pulse_total_4.profile_err = Pulse_total_4.profile()
Pulse_total_4.ph, Pulse_total_4.profilenorm, Pulse_total_4.profile_err = Pulse_total_4.profile_norm()
ph2, ffit_4,A_4,F_4,Sigma_4 = Pulse_total_4.adjusment(5)


#5
Pulse_total_5=pulse_profile()
Pulse_total_5.times=E5
Pulse_total_5.phase_profile=phase_profile
Pulse_total_5.nbin=nbin
Pulse_total_5.T_star_stop=T_star_stop
Pulse_total_5.ph, Pulse_total_5.profile, Pulse_total_5.profile_err = Pulse_total_5.profile()
Pulse_total_5.ph, Pulse_total_5.profilenorm, Pulse_total_5.profile_err = Pulse_total_5.profile_norm()
ph2, ffit_5,A_5,F_5,Sigma_5 = Pulse_total_5.adjusment(5)
#6
Pulse_total_6=pulse_profile()
Pulse_total_6.times=E6
Pulse_total_6.phase_profile=phase_profile
Pulse_total_6.nbin=nbin
Pulse_total_6.T_star_stop=T_star_stop
Pulse_total_6.ph, Pulse_total_6.profile, Pulse_total_6.profile_err = Pulse_total_6.profile()
Pulse_total_6.ph, Pulse_total_6.profilenorm, Pulse_total_6.profile_err = Pulse_total_6.profile_norm()
ph2, ffit_6,A_6,F_6,Sigma_6 = Pulse_total_6.adjusment(5)
#7
Pulse_total_7=pulse_profile()
Pulse_total_7.times=E7
Pulse_total_7.phase_profile=phase_profile
Pulse_total_7.nbin=nbin
Pulse_total_7.T_star_stop=T_star_stop
Pulse_total_7.ph, Pulse_total_7.profile, Pulse_total_7.profile_err = Pulse_total_7.profile()
Pulse_total_7.ph, Pulse_total_7.profilenorm, Pulse_total_7.profile_err = Pulse_total_7.profile_norm()
ph2, ffit_7,A_7,F_7,Sigma_7 = Pulse_total_7.adjusment(5)
#8
Pulse_total_8=pulse_profile()
Pulse_total_8.times=E8
Pulse_total_8.phase_profile=phase_profile
Pulse_total_8.nbin=nbin
Pulse_total_8.T_star_stop=T_star_stop
Pulse_total_8.ph, Pulse_total_8.profile, Pulse_total_8.profile_err = Pulse_total_8.profile()
Pulse_total_8.ph, Pulse_total_8.profilenorm, Pulse_total_8.profile_err = Pulse_total_8.profile_norm()
ph2, ffit_8,A_8,F_8,Sigma_8 = Pulse_total_8.adjusment(5)


plt.figure(figsize=(22.0,7.0))
plt.subplots_adjust(left=0.06, bottom=0.05, right=0.94, top=None, wspace=0.3, hspace=None)
plt.subplot(2,4,1)
plt.suptitle('Source:'+source+'  \n NuSTAR observations \n Pulse period '+str(period)+'s',fontsize=12)

#plt.step(Pulse_total_1.ph,Pulse_total_1.profilenorm,where='mid',color='b',label='3-4.5 KeV')
#plt.errorbar(Pulse_total_1.ph,Pulse_total_1.profilenorm,yerr=Pulse_total_1.profile_err,fmt='bo',markersize=0.5)

plt.plot(ph2,ffit_1(ph2), 'k',label='3-4.5 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_1[0](ph2), 'b--')
plt.plot(ph2,ffit_1[1](ph2), 'r--')
plt.plot(ph2,ffit_1[2](ph2), 'g--')
plt.plot(ph2,ffit_1[3](ph2), 'm--')
plt.plot(ph2,ffit_1[4](ph2), 'y--')
plt.ylim(-15,15)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,2)
#plt.step(Pulse_total_2.ph,Pulse_total_2.profilenorm,where='mid',color='g',label='4.5-6 KeV')
#plt.errorbar(Pulse_total_2.ph,Pulse_total_2.profilenorm,yerr=Pulse_total_2.profile_err,fmt='go',markersize=0.5)

plt.plot(ph2,ffit_2(ph2), 'k',label='4.5-6 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_2[0](ph2), 'b--')
plt.plot(ph2,ffit_2[1](ph2), 'r--')
plt.plot(ph2,ffit_2[2](ph2), 'g--')
plt.plot(ph2,ffit_2[3](ph2), 'm--')
plt.plot(ph2,ffit_2[4](ph2), 'y--')
plt.ylim(-15,15)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,3)
#plt.step(Pulse_total_3.ph,Pulse_total_3.profilenorm,where='mid',color='r',label='6-9 KeV')
#plt.errorbar(Pulse_total_3.ph,Pulse_total_3.profilenorm,yerr=Pulse_total_3.profile_err,fmt='ro',markersize=0.5)

plt.ylim(-15,15)
plt.plot(ph2,ffit_3(ph2), 'k',label='6-9 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_3[0](ph2), 'b--')
plt.plot(ph2,ffit_3[1](ph2), 'r--')
plt.plot(ph2,ffit_3[2](ph2), 'g--')
plt.plot(ph2,ffit_3[3](ph2), 'm--')
plt.plot(ph2,ffit_3[4](ph2), 'y--')
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,4)
#plt.step(Pulse_total_4.ph,Pulse_total_4.profilenorm,where='mid',color='m',label='9-12 KeV')
#plt.errorbar(Pulse_total_4.ph,Pulse_total_4.profilenorm,yerr=Pulse_total_4.profile_err,fmt='mo',markersize=0.5)
plt.plot(ph2,ffit_4(ph2), 'k',label='9-12 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_4[0](ph2), 'b--')
plt.plot(ph2,ffit_4[1](ph2), 'r--')
plt.plot(ph2,ffit_4[2](ph2), 'g--')
plt.plot(ph2,ffit_4[3](ph2), 'm--')
plt.plot(ph2,ffit_4[4](ph2), 'y--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,5)
#plt.step(Pulse_total_5.ph,Pulse_total_5.profilenorm,where='mid',color='b',label='12-15 KeV')
#plt.errorbar(Pulse_total_5.ph,Pulse_total_5.profilenorm,yerr=Pulse_total_5.profile_err,fmt='bo',markersize=0.5)
plt.plot(ph2,ffit_5(ph2), 'k',label='12-15 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_5[0](ph2), 'b--')
plt.plot(ph2,ffit_5[1](ph2), 'r--')
plt.plot(ph2,ffit_5[2](ph2), 'g--')
plt.plot(ph2,ffit_5[3](ph2), 'm--')
plt.plot(ph2,ffit_5[4](ph2), 'y--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,6)
#plt.step(Pulse_total_6.ph,Pulse_total_6.profilenorm,where='mid',color='g',label='15-18 KeV')
#plt.errorbar(Pulse_total_6.ph,Pulse_total_6.profilenorm,yerr=Pulse_total_6.profile_err,fmt='go',markersize=0.5)
plt.plot(ph2,ffit_6(ph2), 'k',label='15-18 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_6[0](ph2), 'b--')
plt.plot(ph2,ffit_6[1](ph2), 'r--')
plt.plot(ph2,ffit_6[2](ph2), 'g--')
plt.plot(ph2,ffit_6[3](ph2), 'm--')
plt.plot(ph2,ffit_6[4](ph2), 'y--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,7)
#plt.step(Pulse_total_7.ph,Pulse_total_7.profilenorm,where='mid',color='r',label='18-27.6 KeV')
#plt.errorbar(Pulse_total_7.ph,Pulse_total_7.profilenorm,yerr=Pulse_total_7.profile_err,fmt='ro',markersize=0.5)
plt.plot(ph2,ffit_7(ph2), 'k',label='18-27.6 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_7[0](ph2), 'b--')
plt.plot(ph2,ffit_7[1](ph2), 'r--')
plt.plot(ph2,ffit_7[2](ph2), 'g--')
plt.plot(ph2,ffit_7[3](ph2), 'm--')
plt.plot(ph2,ffit_7[4](ph2), 'y--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.subplot(2,4,8)
#plt.step(Pulse_total_8.ph,Pulse_total_8.profilenorm,where='mid',color='m',label='27.6-78 KeV')
#plt.errorbar(Pulse_total_8.ph,Pulse_total_8.profilenorm,yerr=Pulse_total_8.profile_err,fmt='mo',markersize=0.5)
plt.plot(ph2,ffit_8(ph2), 'k',label='27.6-78 KeV')#[3,4.5,6,9,12,15,18,27.6,78]
plt.plot(ph2,ffit_8[0](ph2), 'b--')
plt.plot(ph2,ffit_8[1](ph2), 'r--')
plt.plot(ph2,ffit_8[2](ph2), 'g--')
plt.plot(ph2,ffit_8[3](ph2), 'm--')
plt.plot(ph2,ffit_8[4](ph2), 'y--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend()
plt.savefig("NuSTAR_pulse_profile_harmonics.pdf")
plt.show()
#saving the parameters

#A_3,F_3,Sigma_3
AD1=np.concatenate((A_1,F_1,Sigma_1),axis=None)
AD2=np.concatenate((A_2,F_2,Sigma_2),axis=None)
AD3=np.concatenate((A_3,F_3,Sigma_3),axis=None)
AD4=np.concatenate((A_4,F_4,Sigma_4),axis=None)
AD5=np.concatenate((A_5,F_5,Sigma_5),axis=None)
AD6=np.concatenate((A_6,F_6,Sigma_6),axis=None)
AD7=np.concatenate((A_7,F_7,Sigma_7),axis=None)
AD8=np.concatenate((A_8,F_8,Sigma_8),axis=None)
print('popop')
print(AD1)
print(AD2)
print(AD3)
print(AD4)
Parameters_nustar=np.concatenate((AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8),axis=None)
Parameters_nustar=np.array(Parameters_nustar)
print(Parameters_nustar)
print(Parameters_nustar)
hdu = fits.PrimaryHDU(Parameters_nustar)
hdul = fits.HDUList([hdu])
hdul.writeto('Parameters_NUSTAR.fits',overwrite=True)

plt.figure(figsize=(22.0,7.0))
plt.subplot(2,4,1)
plt.suptitle('Source:'+source+'  \n NuSTAR observations \n Pulse period '+str(period)+'s',fontsize=12)
plt.subplots_adjust(left=0.06, bottom=0.05, right=0.94, top=None, wspace=0.3, hspace=None)

plt.step(Pulse_total_1.ph,Pulse_total_1.profilenorm,where='mid',color='b',label='3-4.5 KeV')
plt.errorbar(Pulse_total_1.ph,Pulse_total_1.profilenorm,yerr=Pulse_total_1.profile_err,fmt='bo',markersize=0.5)

plt.plot(ph2,ffit_1(ph2), 'b--')
plt.ylim(-15,15)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,2)
plt.step(Pulse_total_2.ph,Pulse_total_2.profilenorm,where='mid',color='g',label='4.5-6 KeV')
plt.errorbar(Pulse_total_2.ph,Pulse_total_2.profilenorm,yerr=Pulse_total_2.profile_err,fmt='go',markersize=0.5)

plt.plot(ph2,ffit_2(ph2), 'g--')
plt.ylim(-15,15)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,3)
plt.step(Pulse_total_3.ph,Pulse_total_3.profilenorm,where='mid',color='r',label='6-9 KeV')
plt.errorbar(Pulse_total_3.ph,Pulse_total_3.profilenorm,yerr=Pulse_total_3.profile_err,fmt='ro',markersize=0.5)

plt.ylim(-15,15)
plt.plot(ph2,ffit_3(ph2), 'r--')
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,4)
plt.step(Pulse_total_4.ph,Pulse_total_4.profilenorm,where='mid',color='m',label='9-12 KeV')
plt.errorbar(Pulse_total_4.ph,Pulse_total_4.profilenorm,yerr=Pulse_total_4.profile_err,fmt='mo',markersize=0.5)
plt.plot(ph2,ffit_4(ph2), 'm--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,5)
plt.step(Pulse_total_5.ph,Pulse_total_5.profilenorm,where='mid',color='b',label='12-15 KeV')
plt.errorbar(Pulse_total_5.ph,Pulse_total_5.profilenorm,yerr=Pulse_total_5.profile_err,fmt='bo',markersize=0.5)
plt.plot(ph2,ffit_5(ph2), 'b--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,6)
plt.step(Pulse_total_6.ph,Pulse_total_6.profilenorm,where='mid',color='g',label='15-18 KeV')
plt.errorbar(Pulse_total_6.ph,Pulse_total_6.profilenorm,yerr=Pulse_total_6.profile_err,fmt='go',markersize=0.5)
plt.plot(ph2,ffit_6(ph2), 'g--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,7)
plt.step(Pulse_total_7.ph,Pulse_total_7.profilenorm,where='mid',color='r',label='18-27.6 KeV')
plt.errorbar(Pulse_total_7.ph,Pulse_total_7.profilenorm,yerr=Pulse_total_7.profile_err,fmt='ro',markersize=0.5)
plt.plot(ph2,ffit_7(ph2), 'r--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)

plt.subplot(2,4,8)
plt.step(Pulse_total_8.ph,Pulse_total_8.profilenorm,where='mid',color='m',label='27.6-78 KeV')
plt.errorbar(Pulse_total_8.ph,Pulse_total_8.profilenorm,yerr=Pulse_total_8.profile_err,fmt='mo',markersize=0.5)
plt.plot(ph2,ffit_8(ph2), 'm--')
plt.ylim(-18,18)
plt.ylabel('Count rate (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
plt.savefig("NuSTAR_pulse_profile.pdf")
plt.show()

