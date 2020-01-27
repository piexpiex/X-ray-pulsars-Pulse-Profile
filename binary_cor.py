from __future__ import print_function, division
from PyAstronomy import pyasl
import numpy as np

def Binary_orbit(time,asini,porb,ecc,omega_d,t0=-1,t90=-1,pporb=0.0,limit=1.0*10**-6,maxiter=20):
	
	#%\function{BinaryCor}
	#%\synopsis{Removes the influence of the doublestar motion for circular or eliptical orbits.}
	#%\usage{Array_Type t = BinaryCor(Array_Type OR Double_Type time), time in MJD}
	#%\qualifiers{
	#%\qualifier{asini}{: Projected semi-major axis [lt-secs], Mandatory}
	#%\qualifier{porb}{: Orbital period at the epoch [days], Mandatory}
	#%\qualifier{eccentricity}{: Eccentricity, (0<=e<1)}
	#%\qualifier{omega}{: Longitude of periastron [degrees], mandatory}
	#%\qualifier{t0}{: epoch for mean longitude of 0 degrees (periastron, MJD)}
	#%\qualifier{t90}{: epoch for mean longitude of 90 degrees (MJD)}
	#%\qualifier{pporb}{: rate of change of the orbital period (s/s) (default 0)}
	#%\qualifier{limit}{: absolute precision of the correction of the  computation (days, default: 1D-6)}
	#%\qualifier{maxiter}{: stop Newton-Raphson iteration after maxiter steps if limit is not reached (default: 20)}

	#%\description
	#%    (transcribed from IDL-programm BinaryCor.pro)
	#%
	#%    For each time, the z-position of the emitting object is computed
	#%    and the time is adjusted accordingly. This is iterated until
	#%    convergence is reached (usually only one iteration is necessary,
	#%    even in high elliptic cases).
	#%
	#%    Follows equations from Hilditch's book and has also been
	#%    checked against fasebin/axBary. All codes give identical results,
	#%    (to better than 1d-7s) as checked by a Monte Carlo search using
	#%    1d7 different orbits.
	#%
	#%    qualifiers t90 and t0 have to be in days and in the same time
	#3%    system as time (e.g. JD or MJD)
	#%
	#%    Circular orbits:
	#%         * if time of lower conjunction Tlow is known, set
	##%           t0=Tlow and omega=0
	#%         * if time of ascending node is known, Tasc, set
	#%           t90=Tasc and omega=0
	#%         * if time of mid eclipse is known, Tecl, set
	#%           t0=Tecl-0.25*porb and omega=0
	#%!%-

	time=np.array(time)  #make sure that time is an array
	if t90==-1 and t0==-1:
		print("error: need t90 or t0 value")
		return

	if t90!=-1 and t0!=-1:
		print("error: Only one of the t90 and t0 arguments is allowed")
		return

	if ecc<0 :
		print("error: eccentricity must be positive!")
		return
	if ecc>=1:
		print("error: Orbit correction has only been implemented for circular and elliptic orbits")
		return
	if ecc==0:
		omega_d=0 #circular orbit
	
	
	if t0==-1:
		t0 = t90+(omega_d-90.)/360. * porb
		
	if maxiter <=0:
		maxiter=20

	asini_d = asini/86400. #86400 segundos en un día
	t= time

	#Corrections for eccentricity 0<=ecc<1
	omega = omega_d * np.pi/180.0
	sinw = np.sin(omega)
	cosw = np.cos(omega)
	sq = ((1.-ecc)*(1.+ecc))**0.5
	cor =np.array([2.*limit]*len(t))
	

	#start with number of iterations = zero
	numiter=0
	ks = pyasl.MarkleyKESolver() #if you have not PyAstronomy delete this command and use the indications show below
	contada=0
	while((abs(np.amax(cor)) > limit) and (numiter < maxiter)):
		tper = (t-t0)/porb
		m = 2*np.pi*(tper*(1.-0.5*pporb*tper))
		m=np.array(m)
		eanom=np.array([1.0]*len(t))
		#eanom = KeplerEquation(m,ecc)  #if you have not PyAstronomy use this command
		#eanom = KeplerEquation1(m,ecc) #if you have not PyAstronomy use this command (better solution)
		for j in range(0,len(t)):       #if you have not PyAstronomy delete this loop
			eanom[j]=ks.getE(m[j], ecc)       #Mean Anomaly
		sin_e = np.sin(eanom)
		cos_e = np.cos(eanom)
		z = asini_d*(sinw*(cos_e-ecc)+sq*cosw*sin_e)
		f = (t-time)+z                               
		df = (sq*cosw*cos_e - sinw*sin_e)*(2*np.pi*asini_d/(porb*(1.0-ecc*cos_e)))
		cor =f/(1.0+df)
		t = t-cor
		numiter=numiter+1
		contada=contada+1
		print(100*contada/20,"%")
		if numiter >= maxiter:
			print("Exceeded maxiter iterations and did not reach convergence");
			break
	return(t)

def KeplerEquation(m,ecc):
	m=np.array(m)
	if ecc<0 :
		print("error: eccentricity must be positive!")
		return
	if ecc>=1:
		print("error: Orbit correction has only been implemented for circular and elliptic orbits")
	for j in range(0,len(m)):
		
		while m[j]>np.pi:
			m[j]=m[j]-2*np.pi
		while m[j]<-np.pi:
			m[j]=m[j]+2*np.pi
	if ecc==0:
		E=m
	aux=4.0*ecc+0.5
	alpha=(1.0-ecc)/aux
	
	Beta=m/(2.0*aux)
	aux=np.sqrt(Beta**2+alpha**3)
	
	z=Beta+aux
	test=np.array([1.0]*len(z))
	for j in range(0,len(m)):
		if z[j]<=0.0:
			z[j]=beta[j]-aux[j]
			
		test[j]=abs(z[j])**(1/3)
	z=test
	for j in range(0,len(m)):
		if z[j]<0.0:
			z[j]=-z[j]
	s0=z-alpha/z
	
	s1=s0-(0.078*s0**5)/(1.0+ecc)
	e0=m+ecc*(3.0*s1-4.0*s1**3)	
	
	se0=np.sin(e0)
	ce0=np.cos(e0)
	
	f  = e0-ecc*se0-m
	f1 = 1.0-ecc*ce0
	f2 = ecc*se0
	f3 = ecc*ce0
	f4 = -f2
	u1 = -f/f1
	u2 = -f/(f1+0.50*f2*u1)
	u3 = -f/(f1+0.50*f2*u2+.16666666666667*f3*u2*u2)
	u4 = -f/(f1+0.50*f2*u3+.16666666666667*f3*u3*u3+.041666666666667*f4*u3**3)
	
	eccanom=e0+u4
	
	for j in range(0,len(m)):
		while eccanom[j]>=2.0*np.pi:
			eccanom[j]=eccanom[j]-2.0*np.pi
		while eccanom[j]<2.0*np.pi:
			eccanom[j]=eccanom[j]+2.0*np.pi
	
			
	return(eccanom)

def KeplerEquation1(m,ecc):
	m=np.array(m)
	if ecc<0 :
		print("error: eccentricity must be positive!")
		return
	if ecc>=1:
		print("error: Orbit correction has only been implemented for circular and elliptic orbits")
	for j in range(0,len(m)):
		
		while m[j]>np.pi:
			m[j]=m[j]-2*np.pi
		while m[j]<-np.pi:
			m[j]=m[j]+2*np.pi
	if ecc==0:
		E=m
	aux=4.0*ecc+0.5
	alpha=(1.0-ecc)/aux
	
	Beta=m/(2.0*aux)
	aux=np.sqrt(Beta**2+alpha**3)
	
	z=Beta+aux
	test=np.array([1.0]*len(z))
	for j in range(0,len(m)):
		if z[j]<=0.0:
			z[j]=beta[j]-aux[j]
			
		test[j]=abs(z[j])**(1/3)
	z=test
	for j in range(0,len(m)):
		if z[j]<0.0:
			z[j]=-z[j]
	s0=z-alpha/z
	
	s1=s0-(0.078*s0**5)/(1.0+ecc)
	e0=m+ecc*(3.0*s1-4.0*s1**3)	
	
	se0=np.sin(e0)
	ce0=np.cos(e0)
	
	f  = e0-ecc*se0-m
	f1 = 1.0-ecc*ce0
	f2 = ecc*se0
	f3 = ecc*ce0
	f4 = -f2
	u1 = -f/f1
	u2 = -f/(f1+0.50*f2*u1)
	u3 = -f/(f1+0.50*f2*u2+.16666666666667*f3*u2*u2)
	u4 = -f/(f1+0.50*f2*u3+.16666666666667*f3*u3*u3+.041666666666667*f4*u3**3)
	
	eccanom=e0+u4
	
	for j in range(0,len(m)):
		while eccanom[j]>=2.0*np.pi:
			eccanom[j]=eccanom[j]-2.0*np.pi
		while eccanom[j]<2.0*np.pi:
			eccanom[j]=eccanom[j]+2.0*np.pi
	##better solution
	CONT=2
	thresh=10**-5
	for j in range(0,len(m)):
		if m[j]<0:
			m[j]=m[j]+2.0*np.pi
	diff=eccanom-np.sin(eccanom)-m
	contador=0
	for j in range(0,len(m)):
		if abs(diff[j])>10**-10:
			contador=contador+1
	for j in range(0,len(m)):
		if abs(diff[j])>10**-10:
			I=diff[j]
			while CONT>1:
				fe=eccanom[j]-ecc*np.sin(eccanom[j])-m[j]
				fs=1.0-ecc*np.cos(eccanom[j])
				oldval=eccanom[j]
				eccanom[j]=oldval-fe/fs
				if abs(oldval-eccanom[j])<thresh :CONT=0
			while eccanom[j]>= np.pi:eccanom[j]=eccanom[j]-2.0*np.pi
			while eccanom[j]< np.pi:eccanom[j]=eccanom[j]+2.0*np.pi
			
	return(eccanom)
