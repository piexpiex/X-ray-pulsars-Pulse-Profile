import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from read_files import *

READ=read_files()

source=READ[2]
period = READ[8]

#leer los parametros
#XMM
try:
	hdulist = fits.open("Parameters_XMM.fits")
	XMM_parameters = hdulist[0].data
	print(len(XMM_parameters))
except:
	XMM_parameters=np.zeros((60))
#NUSTAR
try:
	hdulist = fits.open("Parameters_NUSTAR.fits")
	NuSTAR_parameters = hdulist[0].data
	print(len(NuSTAR_parameters))
except:
	NuSTAR_parameters=np.zeros((80))
Parameters=np.append(NuSTAR_parameters,XMM_parameters)


A1=[]
A2=[]
A3=[]
A4=[]
A5=[]
F1=[]
F2=[]
F3=[]
F4=[]
F5=[]
DA1=[]
DA2=[]
DA3=[]
DA4=[]
DA5=[]
DF1=[]
DF2=[]
DF3=[]
DF4=[]
DF5=[]
E=[]

#energias

#E=[3,4.5,4.5,6,6,9,9,12,12,15,15,18,18,27.6,27.6,78,0.5,2,2,3.5,3.5,4.5,4.5,6,6,9,9,12]
E=[3,4.5,4.5,6,6,9,9,12,12,15,15,18,18,27.6,27.6,78,0.5,2,2,3,3,4.5,4.5,6,6,9,9,12]

rangos=7
for j in range(int(len(Parameters)/20)):
	A1.append(Parameters[j*20])
	A2.append(Parameters[j*20+1])
	A3.append(Parameters[j*20+2])
	A4.append(Parameters[j*20+3])
	A5.append(Parameters[j*20+4])
	F1.append(Parameters[j*20+5])
	F2.append(Parameters[j*20+6])
	F3.append(Parameters[j*20+7])
	F4.append(Parameters[j*20+8])
	F5.append(Parameters[j*20+9])
	DA1.append(Parameters[j*20+10])
	DA2.append(Parameters[j*20+11])
	DA3.append(Parameters[j*20+12])
	DA4.append(Parameters[j*20+13])
	DA5.append(Parameters[j*20+14])
	DF1.append(Parameters[j*20+15])
	DF2.append(Parameters[j*20+16])
	DF3.append(Parameters[j*20+17])
	DF4.append(Parameters[j*20+18])
	DF5.append(Parameters[j*20+19])

print(A1)
print(A2)
print(A3)
print(A4)
print(A5)

plt.figure(figsize=(22.0,7.0))
plt.subplot(3,2,1)
plt.suptitle('Source:'+source+'  \n Sinusoids amplitudes \n Pulse period '+str(period)+'s',fontsize=12)
plt.subplots_adjust(left=0.06, bottom=0.05, right=0.94, top=None, wspace=None, hspace=0.35)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[A1[j],A1[j]], 'k',label='A1')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A1[j],yerr=DA1[j]**0.5,fmt='ko',markersize=0.5)
	elif j<8:
		plt.plot(E[2*j:2*j+2],[A1[j],A1[j]], 'k')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A1[j],yerr=DA1[j]**0.5,fmt='ko',markersize=0.5)
	else:
		print(E[2*j:2*j+2],[A1[j],A1[j]])
		plt.plot(E[2*j:2*j+2],[A1[j],A1[j]], 'b')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A1[j],yerr=DA1[j]**0.5,fmt='bo',markersize=0.5)
plt.xscale('log')
plt.ylabel('Amplitude')
plt.xlabel("Energy (KeV)")
plt.legend(loc=0)
plt.ylim(0,max(A1)+1)
plt.subplot(3,2,2)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[A2[j],A2[j]], 'k',label='A2')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A2[j],yerr=DA2[j]**0.5,fmt='ko',markersize=0.5)
	elif j<8:
		plt.plot(E[2*j:2*j+2],[A2[j],A2[j]], 'k')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A2[j],yerr=DA2[j]**0.5,fmt='ko',markersize=0.5)
	else:
		plt.plot(E[2*j:2*j+2],[A2[j],A2[j]], 'b')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A2[j],yerr=DA2[j]**0.5,fmt='bo',markersize=0.5)
plt.xscale('log')
plt.ylabel('Amplitude')
plt.xlabel("Energy (KeV)")
plt.legend(loc=2)
plt.ylim(0,max(A2)+1)
plt.subplot(3,2,3)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[A3[j],A3[j]], 'k',label='A3')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A3[j],yerr=DA3[j]**0.5,fmt='ko',markersize=0.5)
	elif j<8:
		plt.plot(E[2*j:2*j+2],[A3[j],A3[j]], 'k')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A3[j],yerr=DA3[j]**0.5,fmt='ko',markersize=0.5)
	else:
		plt.plot(E[2*j:2*j+2],[A3[j],A3[j]], 'b')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A3[j],yerr=DA3[j]**0.5,fmt='bo',markersize=0.5)
plt.xscale('log')
plt.ylabel('Amplitude')
plt.xlabel("Energy (KeV)")
plt.legend(loc=0)
plt.ylim(0,max(A3)+1)
plt.subplot(3,2,4)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[A4[j],A4[j]], 'k',label='A4')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A4[j],yerr=DA4[j]**0.5,fmt='ko',markersize=0.5)
	elif j<8:
		plt.plot(E[2*j:2*j+2],[A4[j],A4[j]], 'k')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A4[j],yerr=DA4[j]**0.5,fmt='ko',markersize=0.5)
	else:
		plt.plot(E[2*j:2*j+2],[A4[j],A4[j]], 'b')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A4[j],yerr=DA4[j]**0.5,fmt='bo',markersize=0.5)
plt.xscale('log')
plt.ylabel('Amplitude')
plt.xlabel("Energy (KeV)")
plt.legend(loc=0)
plt.ylim(0,max(A4)+1)
plt.subplot(3,2,5)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[A5[j],A5[j]], 'k',label='A5')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A5[j],yerr=DA5[j]**0.5,fmt='ko',markersize=0.5)
	elif j<8:
		plt.plot(E[2*j:2*j+2],[A5[j],A5[j]], 'k')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A5[j],yerr=DA5[j]**0.5,fmt='ko',markersize=0.5)
	else:
		plt.plot(E[2*j:2*j+2],[A5[j],A5[j]], 'b')
		plt.errorbar((E[2*j+1]+E[2*j])/2,A5[j],yerr=DA5[j]**0.5,fmt='bo',markersize=0.5)
plt.xscale('log')
plt.ylabel('Amplitude')
plt.xlabel("Energy (KeV)")
plt.legend()
plt.ylim(0,max(A5)+1)
plt.savefig("final.pdf")
plt.show()
#################
#####phase#######
#################

plt.subplot(3,2,1)
plt.suptitle('Source:'+source+'  \n Sinusoids amplitudes \n Pulse period '+str(period)+'s',fontsize=12)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[F1[j],F1[j]], 'k',label='F1')
	elif j<4:
		plt.plot(E[2*j:2*j+2],[F1[j],F1[j]], 'k')
	else:
		plt.plot(E[2*j:2*j+2],[F1[j],F1[j]], 'b')
plt.xscale('log')
plt.ylim(0,1)
plt.legend(loc=2)
plt.subplot(3,2,2)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[F2[j],F2[j]], 'k',label='F2')
	elif j<4:
		plt.plot(E[2*j:2*j+2],[F2[j],F2[j]], 'k')
	else:
		plt.plot(E[2*j:2*j+2],[F2[j],F2[j]], 'b')
plt.xscale('log')
plt.ylim(0,1)
plt.legend(loc=2)
plt.subplot(3,2,3)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[F3[j],F3[j]], 'k',label='F3')
	elif j<4:
		plt.plot(E[2*j:2*j+2],[F3[j],F3[j]], 'k')
	else:
		plt.plot(E[2*j:2*j+2],[F3[j],F3[j]], 'b')
plt.xscale('log')
plt.ylim(0,1)
plt.legend(loc=2)
plt.subplot(3,2,4)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[F4[j],F4[j]], 'k',label='F4')
	elif j<4:
		plt.plot(E[2*j:2*j+2],[F4[j],F4[j]], 'k')
	else:
		plt.plot(E[2*j:2*j+2],[F4[j],F4[j]], 'b')
plt.xscale('log')
plt.ylim(0,1)
plt.legend(loc=2)
plt.subplot(3,2,5)
for j in range(len(A1)):
	if j==0:
		plt.plot(E[2*j:2*j+2],[F5[j],F5[j]], 'k',label='F5')
	elif j<4:
		plt.plot(E[2*j:2*j+2],[F5[j],F5[j]], 'k')
	else:
		plt.plot(E[2*j:2*j+2],[F5[j],F5[j]], 'b')
plt.xscale('log')
plt.ylim(0,1)
plt.legend(loc=2)
plt.show()




