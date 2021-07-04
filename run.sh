: '
#############################
### Files and source name ###
#############################
NuSTAR_file=[ ] # name of the files with NuSTAR data
XMM_file= [ ] # name of the files with XMM-Newton data
source= X-ray pulsar name #source name
#################
### Ephemeris ###
#################
asini= xx   # [It-sec]
Porb= xx   # [days]
ecc= xx      # eccentricity
omega_d= xx  # [degrees]
T0= xx  # [MJD]

pulse_frequency_NuSTAR= 0  # pulsar frequency found in NuSTAR observations (0 to find it or insert the value if known) 
pulse_frequency_XMM= 0 # pulsar frequency found in XMM observations (0 to find it or insert the value if known)
##########################
### working parameters ###
##########################
energy_ranges=[xx,xx,xx] #energy ranges[KeV]
nbin = 40 # number of pulse profile bins
nsinusoids=5 # number of sinusoids for the Fourier series
overwrite= y # (Y/N) [default yes] #overwrite the actual fits files
Z_2_check= n # (Y/N) [default no] #checking the pulse frequency with Z^2 statistic (slows down the analysis)
period_ranges=[xx,xx] # range to search the pulse period [seconds]
period_bins=1000 # number of period bins to test
XMM_Tstart= xx # Initial time of XMM-Newton observation
NuSTAR_Tstart= xx # Initial time of NuSTAR observation
##############################
### running python scripts ###
##############################
'

python python_functions/XMM_data_selector.py
python python_functions/NuSTAR_data_selector.py
python python_functions/XMM_pulse_profile.py
python python_functions/NuSTAR_pulse_profile.py
python python_functions/energy_evolution.py
python python_functions/MERGER.py
