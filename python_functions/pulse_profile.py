import numpy as np
from astropy.modeling import fitting,models
from stingray.pulse.pulsar import fold_events

class pulse_profile():
	times=''
	pulse_frequency=''
	nbin=''
	T_star_stop=''
	ph=''
	profile='' 
	profilenorm='' 
	profile_err = ''
	def profile(self):
		ph, profile, profile_err = fold_events(self.times, self.pulse_frequency, nbin=self.nbin,gtis=self.T_star_stop)
		return(ph, profile, profile_err)
	def profile_norm(self):
		ph=self.ph
		profile=self.profile
		profile_err=self.profile_err
		#normalization of the pulse profile
		mean_profile=0
		for j in range(len(profile)):
			mean_profile=mean_profile+profile[j]
		mean_profile=mean_profile/self.nbin
		std_profile=0
		for j in range(len(profile)):
			std_profile=std_profile+(profile[j]-mean_profile)**2
		std_profile=std_profile**0.5/(self.nbin-1)
		profilenorm=(profile-mean_profile)/std_profile
		profile_err=profile_err/std_profile
		#Duplicating pulse profile
		ph=np.append(ph,ph+1)
		profilenorm=np.append(profilenorm,profilenorm)
		profile_err=np.append(profile_err,profile_err)
		return(ph, profilenorm, profile_err)
		
		
	def adjusment(self,n):
		A=np.zeros(n) #amplitudes
		F=np.zeros(n) #phase origins
		fitter = fitting.LevMarLSQFitter()
		finit = models.Sine1D(frequency=1)
		for j in range(n-1):
			finit=finit +models.Sine1D(frequency=j+2)
		for j in range(n):
			finit[j].frequency.fixed=True
			
		ph2=np.linspace(0,2,200)
		ffit = fitter(finit, self.ph,self.profilenorm, weights = 1/self.profile_err)
		
		for j in range(n):
			A[j]=ffit[j].amplitude.value
			F[j]=ffit[j].phase.value
		for j in range(n):
			while F[j]>1.0:
				F[j]=F[j]-1.0
			while F[j]<0.0:
				F[j]=F[j]+1.0
			if A[j]<0:
				A[j]=-A[j]
				if F[j]>0.5:
					F[j]=F[j]-0.5
				else:
					F[j]=F[j]+0.5
		#Error estimate
		S=0 #chi square
		C_jk= np.zeros((2*n,2*n))
		Sigma=np.zeros(2*n) #frequencies uncertainty
		
		for j in range(self.nbin):
			S=S+(self.profilenorm[j]-ffit(self.ph[j]))**2/self.profile_err[j]**2
		
		for j in range(2*n):
			for k in range(2*n):
				for h in range(self.nbin):
					if j<n and k<n:
						C_jk[j][k]=C_jk[j][k]+ np.sin(2*np.pi*(j+1)*self.ph[h]+F[j])*np.sin(2*np.pi*(k+1)*self.ph[h]+F[k])/self.profile_err[h]**2
					elif j<n and k>n-1:
						C_jk[j][k]=C_jk[j][k]+ np.sin(2*np.pi*(j+1)*self.ph[h]+F[j])*A[k-n]*np.cos(2*np.pi*(k-n+1)*self.ph[h]+F[k-n])/self.profile_err[h]**2
					elif j>n-1 and k<n:
						C_jk[j][k]=C_jk[j][k]+ A[j-n]*np.cos(2*np.pi*(j-n+1)*self.ph[h]+F[j-n])*np.sin(2*np.pi*(k+1)*self.ph[h]+F[k])/self.profile_err[h]**2
					else:
						C_jk[j][k]=C_jk[j][k]+ A[j-n]*np.cos(2*np.pi*(j-n+1)*self.ph[h]+F[j-n])*A[k-n]*np.cos(2*np.pi*(k-n+1)*self.ph[h]+F[k-n])/self.profile_err[h]**2
		
		for j in range(2*n):
			Sigma[j]=(S*np.linalg.inv(C_jk)[j][j]/(self.nbin-n))**0.5
		Sigma_A=Sigma[0:n]
		Sigma_F=Sigma[n:2*n]
		return(ph2,ffit,A,F,Sigma_A,Sigma_F)
