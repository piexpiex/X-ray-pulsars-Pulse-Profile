#############################
### Files and source name ###
#############################

NuSTAR_file=[nu30202004002A01_cl_bary_src.evt,nu30202004002B01_cl_bary_src.evt] # name of the files with the data
XMM_file= 0784570201_evt_bary_src.fits
source= SMC X-1

#################
### Ephemeris ###
#################

asini=53.4876   # [It-sec]
Porb= 3.891923160   # [days]
ecc=0.00089      # eccentricity
omega_d=166   # [degrees]
T0= 52846.688810  # [MJD]

period =0.69965 # pulsar aproximate spin period

pulse_frequency_NuSTAR=1.4280694984937028 # pulsar frequency found in NuSTAR observations(0 to find it or insert the value if known) 
pulse_frequency_XMM=1.4282544464305302 # pulsar frequency found in XMM observations (0 to find it or insert the value if known)

##########################
### working parameters ###
##########################

bin_time = 0.01

nbin = 40

nsinusoids=5

overwrite= n # (Y/N) [default yes] #overwrite the actual fits files

Z_2_check= n # (Y/N) [default no] #checking the pulse frequency with Z^2 statistic (slows down the analysis)

######################
### running python ###
######################

