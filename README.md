# X-ray-pulsars-Pulse-Profile

Program composed by Python and Bash scripts used to obtain pulse profiles of a X-ray pulsar with NuSTAR and XMM-Newton data.

These programs include:

-Select data from Nustar or XMM-Newton fits files

-Search pulse frequency of the pulsar

-Represent the pulse profiles in differents energy ranges

-Adjusment of pulse profiles by five sinusoids (the number of sinusoids is a initial parameter that is possible to change)

-Represent the dependence of the sinusoids parameters with the energy

## Requierements

It is necesary to install *[Stingray][1]* and *[Astropy][2]* Python packages.

[1]: https://stingray.readthedocs.io/en/latest/
[2]: https://www.astropy.org/

## Bash script

run.sh --> you can write the information of the source (ephemeris, name, etc) and the working parameters, later you can run this file in terminal with the command sh run.sh and execute the program.

## Python scripts

binary_cor.py --> Corrects the time values due to the binary orbit delay.

pulse_profile.py --> Define the pulse_profiles class to calculate each pulse profile and perform the fit by sinusoids.

read_files.py --> functions to read and write the information of run.sh.

MERGER.py --> Merge the pdf plots in only one pdf.

###  XMM_data_selector.py / NuSTAR_data_selector.py

Selects the time of the data fits and return fits files with the corrected time values using binary_cor.py.

XMM_data_selector also selects the XMM-Newton PI channel values and return a fits file with this for a faster use by another programs.

This programs plot the Z<sup>2</sup> and the EF statistics variation with a range of frequencies to search the correct frequency.

If the highest best frequency found has a small statistic values, the program ask if you want to search a better frequency by doubling the frequency range (this aspect depends on how exact is the inserted period value).

###  XMM_pulse_profile.py / NuSTAR_pulse_profile.py

Chooses photons by energy ranges and plot their pulse profiles using pulse_profile.py.

Makes an adjustment of the pulse profiles with sinusoids and return a fits file with their parameters and uncertainties.

##  energy_evolution.py

Represents the variation of the fourier parameters with the energy.
