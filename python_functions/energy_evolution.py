import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from read_files import *

READ=read_files()

source=READ[2]
period = READ[8] 

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

Energy_ranges = READ[11]

Energy_ranges=np.array(Energy_ranges)
Energy_ranges=Energy_ranges.astype(np.float)
Energy_ranges=Energy_ranges[np.where((Energy_ranges>=0.5) & (Energy_ranges<=78))]

nsinusoids=READ[13]

#############################
### read parameters files ###
#############################

#XMM-Newton

key_XMM=1
try:
	hdulist = fits.open("fits_folder/"+add_space(source)+"_Parameters_XMM.fits")
	XMM_parameters = hdulist['PARAMETERS'].data
except:
	print('no XMM-Newton data')
	key_XMM=0

#NuSTAR

key_nustar=1
try:
	hdulist = fits.open("fits_folder/"+add_space(source)+"_Parameters_NUSTAR.fits")
	NuSTAR_parameters = hdulist['PARAMETERS'].data
except:
	print('no NuSTAR data')
	key_nustar=0
	

########################################
### Plot parameters energy evolution ###
########################################

# Amplitudes

plt.subplots_adjust(left=0.06, bottom=0.05, right=0.94, top=None, wspace=None, hspace=0.35)
plt.figure(figsize=(22.0,7.0))

for k in range(nsinusoids):
	plt.subplot(round((nsinusoids+0.1)/2),2,k+1)
	#plt.suptitle('Source:'+source+'  \n Sinusoids amplitudes \n Pulse period '+str(period)+'s',fontsize=12)
	max_value=0
	#NuSTAR
	if key_nustar==1:
		for j in range(len(NuSTAR_parameters)):
			if j==0:
				plt.plot([NuSTAR_parameters[j][1],NuSTAR_parameters[j][2]],[NuSTAR_parameters[j][3+k*4],NuSTAR_parameters[j][3+k*4]], 'k',label='A'+str(k+1))
				plt.errorbar((NuSTAR_parameters[j][1]+NuSTAR_parameters[j][2])/2,NuSTAR_parameters[j][3+k*4],yerr=NuSTAR_parameters[j][5+k*4]**0.5,fmt='ko',markersize=0.5)
			else:
				plt.plot([NuSTAR_parameters[j][1],NuSTAR_parameters[j][2]],[NuSTAR_parameters[j][3+k*4],NuSTAR_parameters[j][3+k*4]], 'k')
				plt.errorbar((NuSTAR_parameters[j][1]+NuSTAR_parameters[j][2])/2,NuSTAR_parameters[j][3+k*4],yerr=NuSTAR_parameters[j][5+k*4]**0.5,fmt='ko',markersize=0.5)
			if NuSTAR_parameters[j][3+k*4]>max_value:
				max_value=NuSTAR_parameters[j][3+k*4]
	#XMM-Newton
	if key_XMM==1:
		for j in range(len(XMM_parameters)):
			if key_nustar==0 and j==0:
				plt.plot([XMM_parameters[j][1],XMM_parameters[j][2]],[XMM_parameters[j][3+k*4],XMM_parameters[j][3+k*4]], 'b',label='A'+str(k+1))
				plt.errorbar((XMM_parameters[j][1]+XMM_parameters[j][2])/2,XMM_parameters[j][3+k*4],yerr=XMM_parameters[j][5+k*4]**0.5,fmt='bo',markersize=0.5)
			else:
				plt.plot([XMM_parameters[j][1],XMM_parameters[j][2]],[XMM_parameters[j][3+k*4],XMM_parameters[j][3+k*4]], 'b')
				plt.errorbar((XMM_parameters[j][1]+XMM_parameters[j][2])/2,XMM_parameters[j][3+k*4],yerr=XMM_parameters[j][5+k*4]**0.5,fmt='bo',markersize=0.5)
			if XMM_parameters[j][3+k*4]>max_value:
				max_value=XMM_parameters[j][3+k*4]
	
	plt.xscale('log')
	plt.ylabel('Amplitude')
	plt.xlabel("Energy (KeV)")
	plt.ylim(0,max_value+1)
	plt.legend()


plt.savefig('figures_folder/'+add_space(source)+'_amplitudes_energy_variation.pdf')
