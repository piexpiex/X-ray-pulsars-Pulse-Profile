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
from astropy.modeling import fitting,models

name_file="src_sd.events"

hdulist = fits.open(name_file)
hdulist.info()
header = hdulist[2].header
print('...')
#print(header)
print('...')

#datos (en principio)
asini=26.33 #[lt-secs]
Porb=1.4084 #[days]
ecc=0.0     
omega_d=0.0 #[degrees]
T0=51110.866#(MJD)

period = 13.5 #s

bin_time = 0.01

nbin = 40

hdulist = fits.open(name_file)
hdulist.info()
header = hdulist[0].header

#lo anterior
ev =EventList()
ev = ev.read(name_file, 'fits')


Tstart=hdulist['STDGTI04'].data['START']
Tstop=hdulist['STDGTI04'].data['STOP']
#Tstart=hdulist['GTI00007'].data['START']
#Tstop=hdulist['GTI00007'].data['STOP']
print('T start',Tstart)
CT=ev.ncounts
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
print('T_star_stop',T_star_stop)
phase_eso=0.074131 #0.074185
print('tiempo total=',Total_exptime)
ph, profile, profile_err = fold_events(times,phase_eso, nbin=nbin,gtis=T_star_stop)#0.00356

obs_length = times[len(times)-1]-times[0]

#la parte de energías #0.0035298

#E1 95-600
A1=np.where((PI<600) & (PI>95))
E1=times[A1]
pE1, profileE, profile_errE = fold_events(E1,phase_eso, nbin=nbin,gtis=T_star_stop)
#plt.plot(ph, profileE)
#E1 600-1200
A2=np.where((PI<1200) & (PI>600))
E2=times[A2]
pE2, profileE2, profile_errE2 = fold_events(E2,phase_eso, nbin=nbin,gtis=T_star_stop)
#plt.plot(ph, profileE2)
#E1 1200-2400
A3=np.where((PI<2400) & (PI>1200))
E3=times[A3]
print(len(E3))
print('koko',len(E3))
pE3, profileE3, profile_errE3 = fold_events(E3,phase_eso, nbin=nbin,gtis=T_star_stop)
#plt.plot(ph, profileE3)


Err=profile_errE3*nbin/Total_exptime
Err1=profile_errE3*nbin/Total_exptime
Err2=profile_errE3*nbin/Total_exptime
Err3=profile_errE3*nbin/Total_exptime
profile=profile*nbin/Total_exptime
profileE=profileE*nbin/Total_exptime
profileE2=profileE2*nbin/Total_exptime
profileE3=profileE3*nbin/Total_exptime

#normalization and adjusment
fitter = fitting.LevMarLSQFitter()

mean_profile=0
for j in range(len(profile)):
	mean_profile=mean_profile+profile[j]
mean_profile=mean_profile/nbin
std_profile=0
for j in range(len(profile)):
	std_profile=std_profile+(profile[j]-mean_profile)**2
std_profile=std_profile**0.5/(nbin-1)
profilenorm=(profile-mean_profile)/std_profile
Err=Err/std_profile
#E1
mean_profile1=0
for j in range(len(profile)):
	mean_profile1=mean_profile1+profileE[j]
mean_profile1=mean_profile1/nbin
std_profile1=0
for j in range(len(profile)):
	std_profile1=std_profile1+(profileE[j]-mean_profile1)**2
std_profile1=std_profile1**0.5/(nbin-1)
profilenormE=(profileE-mean_profile1)/std_profile1
Err1=Err1/std_profile1
#E2
mean_profile2=0
for j in range(len(profile)):
	mean_profile2=mean_profile2+profileE2[j]
mean_profile2=mean_profile2/nbin
std_profile2=0
for j in range(len(profile)):
	std_profile2=std_profile2+(profileE2[j]-mean_profile2)**2
std_profile2=std_profile2**0.5/(nbin-1)
profilenormE2=(profileE2-mean_profile2)/std_profile2
Err2=Err2/std_profile2
#E3
mean_profile3=0
for j in range(len(profile)):
	mean_profile3=mean_profile3+profileE3[j]
mean_profile3=mean_profile3/nbin
std_profile3=0
for j in range(len(profile)):
	std_profile3=std_profile3+(profileE3[j]-mean_profile3)**2
std_profile3=std_profile3**0.5/(nbin-1)
profilenormE3=(profileE3-mean_profile3)/std_profile3
Err3=Err3/std_profile3

ph=np.append(ph,ph+1)
ph2=np.linspace(0,2,200)
profilenorm=np.append(profilenorm,profilenorm)
profilenormE=np.append(profilenormE,profilenormE)
profilenormE2=np.append(profilenormE2,profilenormE2)
profilenormE3=np.append(profilenormE3,profilenormE3)
Err=np.append(Err,Err)
Err1=np.append(Err1,Err1)
Err2=np.append(Err2,Err2)
Err3=np.append(Err3,Err3)
#plt.text(4, 9, "Source:1E1145.1-6141  \n NuSTAR observations \n Pulse period 295s", size=10,bbox=dict(boxstyle="round",fc="w", ec="k"))

Parameters_nustar=[]
plt.subplot(2,2,1)

#Error estimate
S=np.zeros(4) #statistical deviation of any adjustment and its real values
Ckk_A= np.zeros((4,5)) #diagonal elements of C matrix (amplitudes)
Ckk_F=np.zeros((4,5)) #diagonal elements of C matrix (initial phase value)
Sigma_A=np.zeros((4,5)) #amplitudes uncertainty
Sigma_F=np.zeros((4,5)) #frequencies uncertainty


finit = models.Const1D(amplitude=0)+models.Sine1D(frequency=1)+models.Sine1D(frequency=2)+models.Sine1D(frequency=3)+models.Sine1D(frequency=4)+models.Sine1D(frequency=5)
finit[0].amplitude.fixed=True
finit[1].frequency.fixed=True
finit[2].frequency.fixed=True
finit[3].frequency.fixed=True
finit[4].frequency.fixed=True
finit[5].frequency.fixed=True

plt.suptitle('Source:LMC X-4  \n XMM-Newton observations \n Pulse period 13.5s',fontsize=12)
plt.step(ph,profilenorm,where='mid',color='k',label='total pulse')
plt.errorbar(ph,profilenorm,yerr=Err,fmt='ko',markersize=0.5)
ffit = fitter(finit, ph,profilenorm, weights = 1/Err1)
print(ffit[1].frequency)
print(ffit[1])
print(ffit[2])
print(ffit[3])
print(ffit[4])
print(ffit[5])
print('---total---')
print(ffit[0].amplitude.value)
print(ffit[1].amplitude.value,ffit[1].phase.value)
print(ffit[2].amplitude.value,ffit[2].phase.value)
print(ffit[3].amplitude.value,ffit[3].phase.value)
print(ffit[4].amplitude.value,ffit[4].phase.value)
print(ffit[5].amplitude.value,ffit[5].phase.value)
Parameters_nustar.append(ffit[0].amplitude.value)
Parameters_nustar.append(ffit[1].amplitude.value)
Parameters_nustar.append(ffit[1].phase.value)
Parameters_nustar.append(ffit[2].amplitude.value)
Parameters_nustar.append(ffit[2].phase.value)
Parameters_nustar.append(ffit[3].amplitude.value)
Parameters_nustar.append(ffit[3].phase.value)
Parameters_nustar.append(ffit[4].amplitude.value)
Parameters_nustar.append(ffit[4].phase.value)
Parameters_nustar.append(ffit[5].amplitude.value)
Parameters_nustar.append(ffit[5].phase.value)
plt.ylim(-15,15)
plt.plot(ph2,ffit(ph2), 'k--')
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
for j in range(nbin):
	S[0]=S[0]+(profilenorm[j]-ffit((2*j+1)/nbin))**2/Err[j]**2
#S[0]=S[0]**0.5/(nbin-1)
for j in range(5):
	for k in range(nbin):
		Ckk_F[0][j]=Ckk_F[0][j]+ffit[j+1].amplitude.value**2*np.cos(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err[k]**2
		Ckk_A[0][j]=Ckk_A[0][j]+np.sin(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err[k]**2
for j in range(5):
	Sigma_A[0][j]=S[0]/(nbin-10)/Ckk_A[0][j]
	Sigma_F[0][j]=S[0]/(nbin-10)/Ckk_F[0][j]
plt.subplot(2,2,2)
plt.step(ph,profilenormE,where='mid',color='b',label='0,5-3 KeV')
plt.errorbar(ph,profilenormE,yerr=Err1,fmt='bo',markersize=0.5)
ffit = fitter(finit, ph,profilenormE, weights = 1/Err1)
print('---0,5-3,0 eV---')
print(ffit[0].amplitude.value)
print(ffit[1].amplitude.value,ffit[1].phase.value)
print(ffit[2].amplitude.value,ffit[2].phase.value)
print(ffit[3].amplitude.value,ffit[3].phase.value)
print(ffit[4].amplitude.value,ffit[4].phase.value)
print(ffit[5].amplitude.value,ffit[5].phase.value)
Parameters_nustar.append(ffit[0].amplitude.value)
Parameters_nustar.append(ffit[1].amplitude.value)
Parameters_nustar.append(ffit[1].phase.value)
Parameters_nustar.append(ffit[2].amplitude.value)
Parameters_nustar.append(ffit[2].phase.value)
Parameters_nustar.append(ffit[3].amplitude.value)
Parameters_nustar.append(ffit[3].phase.value)
Parameters_nustar.append(ffit[4].amplitude.value)
Parameters_nustar.append(ffit[4].phase.value)
Parameters_nustar.append(ffit[5].amplitude.value)
Parameters_nustar.append(ffit[5].phase.value)
plt.plot(ph2,ffit(ph2), 'b--')
plt.ylim(-15,15)
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
for j in range(nbin):
	S[1]=S[1]+(profilenormE[j]-ffit((2*j+1)/nbin))**2/Err1[j]**2
#S[1]=S[1]**0.5/(nbin-1)
for j in range(5):
	for k in range(nbin):
		Ckk_F[1][j]=Ckk_F[1][j]+ffit[j+1].amplitude.value**2*np.cos(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err1[k]**2
		Ckk_A[1][j]=Ckk_A[1][j]+np.sin(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err1[k]**2
for j in range(5):
	Sigma_A[1][j]=S[1]/(nbin-10)/Ckk_A[1][j]
	Sigma_F[1][j]=S[1]/(nbin-10)/Ckk_F[1][j]
plt.subplot(2,2,3)
plt.step(ph,profilenormE2,where='mid',color='g',label='3-6 KeV')
plt.errorbar(ph,profilenormE2,yerr=Err2,fmt='go',markersize=0.5)
ffit = fitter(finit, ph,profilenormE2, weights = 1/Err1)
print('---3,0-6,0 eV---')
print(ffit[0].amplitude.value)
print(ffit[1].amplitude.value,ffit[1].phase.value)
print(ffit[2].amplitude.value,ffit[2].phase.value)
print(ffit[3].amplitude.value,ffit[3].phase.value)
print(ffit[4].amplitude.value,ffit[4].phase.value)
print(ffit[5].amplitude.value,ffit[5].phase.value)
Parameters_nustar.append(ffit[0].amplitude.value)
Parameters_nustar.append(ffit[1].amplitude.value)
Parameters_nustar.append(ffit[1].phase.value)
Parameters_nustar.append(ffit[2].amplitude.value)
Parameters_nustar.append(ffit[2].phase.value)
Parameters_nustar.append(ffit[3].amplitude.value)
Parameters_nustar.append(ffit[3].phase.value)
Parameters_nustar.append(ffit[4].amplitude.value)
Parameters_nustar.append(ffit[4].phase.value)
Parameters_nustar.append(ffit[5].amplitude.value)
Parameters_nustar.append(ffit[5].phase.value)
plt.plot(ph2,ffit(ph2), 'g--')
plt.ylim(-15,15)
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
for j in range(nbin):
	S[2]=S[2]+(profilenormE[j]-ffit((2*j+1)/nbin))**2/Err2[j]**2
#S[2]=S[2]**0.5/(nbin-1)
for j in range(5):
	for k in range(nbin):
		Ckk_F[2][j]=Ckk_F[2][j]+ffit[j+1].amplitude.value**2*np.cos(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err2[k]**2
		Ckk_A[2][j]=Ckk_A[2][j]+np.sin(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err2[k]**2
for j in range(5):
	Sigma_A[2][j]=S[2]/(nbin-10)/Ckk_A[2][j]
	Sigma_F[2][j]=S[2]/(nbin-10)/Ckk_F[2][j]
plt.subplot(2,2,4)
plt.step(ph,profilenormE3,where='mid',color='r',label='6-12 KeV')
plt.errorbar(ph,profilenormE3,yerr=Err3,fmt='ro',markersize=0.5)
ffit = fitter(finit,ph,profilenormE3, weights = 1/Err1)
print('---6,0-12 eV---')
print(ffit[0].amplitude.value)
print(ffit[1].amplitude.value,ffit[1].phase.value)
print(ffit[2].amplitude.value,ffit[2].phase.value)
print(ffit[3].amplitude.value,ffit[3].phase.value)
print(ffit[4].amplitude.value,ffit[4].phase.value)
print(ffit[5].amplitude.value,ffit[5].phase.value)
Parameters_nustar.append(ffit[0].amplitude.value)
Parameters_nustar.append(ffit[1].amplitude.value)
Parameters_nustar.append(ffit[1].phase.value)
Parameters_nustar.append(ffit[2].amplitude.value)
Parameters_nustar.append(ffit[2].phase.value)
Parameters_nustar.append(ffit[3].amplitude.value)
Parameters_nustar.append(ffit[3].phase.value)
Parameters_nustar.append(ffit[4].amplitude.value)
Parameters_nustar.append(ffit[4].phase.value)
Parameters_nustar.append(ffit[5].amplitude.value)
Parameters_nustar.append(ffit[5].phase.value)
plt.plot(ph2,ffit(ph2), 'r--')
plt.ylim(-15,15)
plt.ylabel('Counts/sec (normalized)')
plt.xlabel("$\phi$")
plt.legend(loc=2)
for j in range(nbin):
	S[3]=S[3]+(profilenormE[j]-ffit((2*j+1)/nbin))**2/Err3[k]**2
#S[3]=S[3]**0.5/(nbin-1)
for j in range(5):
	for k in range(nbin):
		Ckk_F[3][j]=Ckk_F[3][j]+ffit[j+1].amplitude.value**2*np.cos(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err3[k]**2
		Ckk_A[3][j]=Ckk_A[3][j]+np.sin(2*np.pi*k*(2*k+1)/nbin +ffit[j+1].phase.value )**2/Err3[k]**2
for j in range(5):
	Sigma_A[3][j]=S[3]/(nbin-10)/Ckk_A[3][j]
	Sigma_F[3][j]=S[3]/(nbin-10)/Ckk_F[3][j]
plt.show()

print('---uncertainties---')
print(' ')
print('---Amplitudes---')
print(Sigma_A)
print(' ')
print('---Frequencies---')
print(Sigma_F)
#guardar los parametros
Parameters_nustar=np.array(Parameters_nustar)
print(Parameters_nustar)
hdu = fits.PrimaryHDU(Parameters_nustar)
hdul = fits.HDUList([hdu])
hdul.writeto('Parameters_XMM.fits')
