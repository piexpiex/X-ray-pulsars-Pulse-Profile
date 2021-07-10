# X-ray-pulsars-Pulse-Profile

Program composed by Python and Bash scripts used to obtain pulse profiles of a X-ray pulsar with NuSTAR and XMM-Newton data.

To use it modify the run.sh parameters and use the command "sh run.sh" with the python_functions folder and data files in fits format in the same directory.

These programs include:

-Select data from Nustar or XMM-Newton (EPIC-pn) fits files

-Data analysis and technical processes

-Search pulse frequency of the pulsar in each observation with epoch folding distribution (and if it were possible the epoch folding distribution would be fitted by a Lorentzian curve in order to obtain the respective uncertainty)

-Represent the pulse profiles in differents energy ranges
  
-Fit the pulse profiles by a number of sinusoids

-Represent the variation of the sinusoids parameters with the energy
  
## Requierements

It is necesary to install *[Stingray][1]* and *[Astropy][2]* Python packages.

[1]: https://stingray.readthedocs.io/en/latest/
[2]: https://www.astropy.org/

If you want to merge the pdfs you have to also have *[PyPDF2][3]* python library to be able to use MERGER.py script.

[3]: https://pythonhosted.org/PyPDF2/

## Bash script

run.sh --> you can write the information of the source (ephemeris, name, etc) and the working parameters, later you can run this file in terminal with the command sh run.sh and execute the program.

-------------------------------
                WARNING
The run.sh file may not work properly on some operating systems like Windows, if this happens try rewriting the file or writing the following commands in Unix terminal:

$ python python_functions/XMM_data_selector.py <br/>
$ python python_functions/NuSTAR_data_selector.py <br/>
$ python python_functions/XMM_pulse_profile.py <br/>
$ python python_functions/NuSTAR_pulse_profile.py <br/>
$ python python_functions/energy_evolution.py <br/>
$ python python_functions/MERGER.py

--------------------------------------

## Python scripts

binary_cor.py --> Corrects the time values due to the binary orbit delay (translated functions of *[IAAT][4]*).

[4]: http://astro.uni-tuebingen.de/software/idl/aitlib/astro/

pulse_profile.py --> Define the pulse_profiles class to calculate each pulse profile and perform the fit by sinusoids.

read_files.py --> functions to read and write the information of run.sh.

MERGER.py --> Merge the pdf plots in only one summary pdf (using PyPDF2).

###  XMM_data_selector.py / NuSTAR_data_selector.py

selects the time of the data fits and return fits files with the corrected time values using binary_cor.py.

XMM_data_selector also selects the XMM-Newton PI channel values and return a fits file with this for a faster use by another programs.

plots the Z<sup>2</sup> and the EF statistics variation with a range of frequencies to search the correct frequency.

If the highest best frequency found has a small statistic values, the program ask if you want to search a better frequency by doubling the frequency range (this aspect depends on how exact is the inserted period value).

###  XMM_pulse_profile.py / NuSTAR_pulse_profile.py

colects photons of the standar good time intervals (STDGTI)

classifies photons by energy ranges and plot their pulse profiles using pulse_profile.py.

makes an adjustment of the pulse profiles with sinusoids and return a fits file with their parameters and uncertainties and other with the pulses and their fits.

###  energy_evolution.py

represents the variation of the fourier parameters with the energy.

## run.sh parameters

run.sh script currently has 17 editable parameters, these are described below:

NuSTAR_file=[NuSTAR_file_1.fits,NuSTAR_file_2.fits,...] --> name of the files with NuSTAR data in fits format write in square brackets (if you do not have NuSTAR data put "[ ]")

XMM_file= [XMM_file_1.fits,XMM_file_2.fits,...] --> name of the files with XMM-Newton data in fits format write in square brackets (if you do not have XMM_Newton data put "[ ]")

source= X-ray pulsar name --> name of the analized source  
asini= xx --> projected semi-major axis in It-sec  
Porb= xx --> orbital period at the epoch in days  
ecc= xx -->  eccentricity  
omega_d= xx --> longitude of periastron in degrees  
T0= xx --> epoch for mean longitude of 0 degrees in MJD system

(to perform the analysis without the binay orbit correction, set asini=0, Porb=0, ecc=0, omega_d=0 and T0=0)

pulse_frequency_NuSTAR= 0 --> pulsar frequency found in NuSTAR observations (0 to find it or insert the value if known) 
pulse_frequency_XMM= 0 --> pulsar frequency found in XMM-Newton observations (0 to find it or insert the value if known)

energy_ranges=[xx,xx,xx] --> energy ranges in KeV (XMM_Newton --> 0.5-12 & NuSTAR --> 3-78) in square brackets

nbin = xx --> number of pulse profile bins

nsinusoids= xx --> number of sinusoids for the Fourier series fit

overwrite= xx --> (Y/N) [default yes] command to overwrite the actual fits files obtained by XMM_data_selector.py/NuSTAR_data_selector.py, if you want to perform the analysis in multiple stages select "N" to do not repeat the first stages, but you can continue to change the energy ranges or another analysis parameter

Z_2_check= xx --> (Y/N) [default no] command to check the pulse frequency with Z^2 statistic (slows down the analysis)

period_ranges=[xx,xx] --> range to search the pulse period

period_bins= xx --> number of period bins to test

XMM_Tstart= xx --> Initial time of XMM-Newton observation  <br/>
NuSTAR_Tstart= xx --> Initial time of NuSTAR observation
  


