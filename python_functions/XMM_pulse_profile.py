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
import os

READ=read_files() # run.sh data
name_file=READ[1] 

if name_file[0]==' ':
	exit()

source=READ[2]

period = READ[8]
XMM_time=READ[18]

if READ[9]!=0 or READ[10]!=0:
	if READ[9]!=0 or READ[10]==0:
		period=1/READ[9]
		period =round(period,2-int(np.log10(period)))
	elif READ[9]==0 or READ[10]!=0:
		period=1/READ[10]
		period =round(period,2-int(np.log10(period)))
	else:
		period=(1/READ[9]+1/READ[10])/2
		period =round(period,2-int(np.log10(period)))
		
else:
	print('no avalaible spin frequency')
	exit()

pulse_frequency=READ[10] 

Energy_ranges = READ[11]

nbin = READ[12]

nsinusoids=READ[13]

try:
	os.mkdir('fits_folder')
	os.mkdir('figures_folder')
	print('making directories')
	print('starting the analysis')
except:
	print('starting the analysis')

####################
### Data lecture ###
####################

Tstart=np.array([])
Tstop=np.array([])
for k in range(len(name_file)):
	hdulist = fits.open(name_file[k])
	header = hdulist[0].header
	for j in range(1,13): # lecture of the standar good time interval is maximal to 12 (equal to the numbre of pn chips)
		if j>9:
			n_stdgti=str(j)
		else:
			n_stdgti='0'+str(j)
		try:
			Tstart=np.append(Tstart,hdulist['STDGTI'+n_stdgti].data['START'])
			Tstop=np.append(Tstop,hdulist['STDGTI'+n_stdgti].data['STOP'])
		except:
			continue



#times & PI
hdulist = fits.open("fits_folder/XMM_times.fits")
data = hdulist['VALUES'].data
times=data.field(0) 
PI =data.field(1)

T_star_stop=[[Tstart[0],Tstop[0]]]
Total_exptime=0
for j in range(len(Tstart)):
	Total_exptime=Total_exptime+Tstop[j]-Tstart[j]
	if j>0:
		T_star_stop.append([Tstart[j],Tstop[j]])
T_star_stop=np.array(T_star_stop)

times=times-XMM_time
T_star_stop=T_star_stop-XMM_time

#####################
### Pulse profile ###
#####################

Energy_ranges=np.array(Energy_ranges)
Energy_ranges=Energy_ranges.astype(np.float)
Energy_ranges=Energy_ranges[np.where((Energy_ranges>=0.5) & (Energy_ranges<=12))]

def PI_calculator(E):
	PI_range=[0.0]*len(E)
	for j in range(len(E)):
		if E[j]<1.0:
			PI_range[j]= -20*(E[j]-1)**2+200*E[j] #little correction for non linearity below 2 KeV
		else:
			PI_range[j]= 0+200*E[j]
	return(np.array(PI_range))

PI_ranges=PI_calculator(Energy_ranges)

pi_photons=[[0.0]]*(len(Energy_ranges)-1) 
time_photons=[[0.0]]*(len(Energy_ranges)-1)

for j in range(len(pi_photons)):
	pi_photons[j]=np.where((PI<PI_ranges[j+1]) & (PI>PI_ranges[j]))
	time_photons[j]=times[pi_photons[j]]

Pulse_profiles=[]
ffit=[]
amplitudes=[]
initial_phases=[]
amplitudes_sigma=[]
initial_phases_sigma=[]

for j in range(len(time_photons)):
	Pulse_profiles.append(pulse_profile())
	Pulse_profiles[j].times=time_photons[j]
	Pulse_profiles[j].pulse_frequency=pulse_frequency
	Pulse_profiles[j].nbin=nbin
	Pulse_profiles[j].T_star_stop=T_star_stop
	
	Pulse_profiles[j].ph, Pulse_profiles[j].profile,Pulse_profiles[j].profile_err = Pulse_profiles[j].profile()
	Pulse_profiles[j].ph, Pulse_profiles[j].profilenorm, Pulse_profiles[j].profile_err = Pulse_profiles[j].profile_norm()	
	ph2, ffit_j,A_j,F_j,Sigma_Aj,Sigma_Fj = Pulse_profiles[j].adjusment(nsinusoids)
	
	ffit.append(ffit_j)
	amplitudes.append(A_j)
	initial_phases.append(F_j)
	amplitudes_sigma.append(Sigma_Aj)
	initial_phases_sigma.append(Sigma_Fj)
	
#saving the parameters



c1 = fits.Column(name='Energy range', array=np.arange(1,len(Energy_ranges)), format='E')
c2 = fits.Column(name='Initial energy (KeV)', array=Energy_ranges[0:len(Energy_ranges)-1], format='E')
c3 = fits.Column(name='Final energy (KeV)', array=Energy_ranges[1:len(Energy_ranges)], format='E')
COLUMNS=[c1,c2,c3]
for j in range(nsinusoids):
	A_j=[]
	F_j=[]
	Sigma_Aj=[]
	Sigma_Fj=[]
	for k in range(len(time_photons)):
		A_j.append(amplitudes[k][j])
		F_j.append(initial_phases[k][j])
		Sigma_Aj.append(amplitudes_sigma[k][j])
		Sigma_Fj.append(initial_phases_sigma[k][j])
	c_A_j=fits.Column(name='A'+str(j+1), array=np.array(A_j), format='D')
	COLUMNS.append(c_A_j)
	c_F_j=fits.Column(name='F'+str(j+1), array=np.array(F_j), format='D')
	COLUMNS.append(c_F_j)
	c_Sigma_Aj=fits.Column(name='Sigma A'+str(j+1), array=np.array(Sigma_Aj), format='D')
	COLUMNS.append(c_Sigma_Aj)
	c_Sigma_Fj=fits.Column(name='Sigma F'+str(j+1), array=np.array(Sigma_Fj), format='D')
	COLUMNS.append(c_Sigma_Fj)


t = fits.BinTableHDU.from_columns(COLUMNS,name='parameters')
t.writeto('fits_folder/Parameters_XMM.fits',overwrite=True)

############################
### Pulse profiles plots ###
############################

plt.figure(figsize=(22.0,7.0))
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.94, top=None, wspace=0.3, hspace=None)
plt.suptitle('Source:'+source+'  \n XMM-Newton observations \n Pulse period '+str(period)+'s',fontsize=12)
for j in range(len(time_photons)):
	if len(time_photons)>2:
		plt.subplot(2,round((len(time_photons)+0.1)/2),j+1)
	else:
		plt.subplot(1,round((len(time_photons)+0.1)),j+1)

	plt.plot(ph2,ffit[j](ph2), 'k',label=str(Energy_ranges[j])+'-'+str(Energy_ranges[j+1])+' KeV')
	for k in range(nsinusoids):
		if k==0:
			plt.plot(ph2,ffit[j][k](ph2), 'b--')
		if k==1:
			plt.plot(ph2,ffit[j][k](ph2), 'r--')
		if k==2:
			plt.plot(ph2,ffit[j][k](ph2), 'g--')
		if k==3:
			plt.plot(ph2,ffit[j][k](ph2), 'm--')
		if k>3:
			plt.plot(ph2,ffit[j][k](ph2), 'y--')
	plt.ylim(min(ffit[j](ph2))*1.2,max(ffit[j](ph2))*1.2)
	plt.ylabel('Count rate (normalized)')
	plt.xlabel("$\phi$")
	plt.legend()

plt.savefig('figures_folder/'+add_space(source)+'_XMM_pulse_profile_harmonics.pdf')

plt.figure(figsize=(22.0,7.0))
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.94, top=None, wspace=0.3, hspace=None)
plt.suptitle('Source:'+source+'  \n XMM-Newton observations \n Pulse period '+str(period)+'s',fontsize=12)

color=['k']*10
for j in range(len(time_photons)):
	if len(time_photons)>2:
		plt.subplot(2,round((len(time_photons)+0.1)/2),j+1)
	else:
		plt.subplot(1,round((len(time_photons)+0.1)),j+1)

	plt.step(Pulse_profiles[j].ph,Pulse_profiles[j].profilenorm,where='mid',color=color[j],label=str(Energy_ranges[j])+'-'+str(Energy_ranges[j+1])+' KeV')
	plt.errorbar(Pulse_profiles[j].ph,Pulse_profiles[j].profilenorm,yerr=Pulse_profiles[j].profile_err,fmt=color[j]+'o',markersize=0.5)

	plt.plot(ph2,ffit[j](ph2), color[j]+'--')

	plt.ylim(min(ffit[j](ph2))*1.2,max(ffit[j](ph2))*1.2)
	plt.ylabel('Count rate (normalized)')
	plt.xlabel("$\phi$")
	plt.legend()

plt.savefig('figures_folder/'+add_space(source)+'_XMM_pulse_profile.pdf')
