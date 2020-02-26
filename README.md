# X-ray-pulsars-Pulse-Profile

Programs used to obtain pulse profiles of X-ray pulsar with NuSTAR and XMM-Newton data.

These programs include:

-Select data from Nustar or XMM-Newton fits files

-Search frequency of the pulsar

-Represent the pulse profiles in differents energy ranges

-Adjusment of pulse profiles by five sinusoids (the number of sinusoids is a initial parameter that is possible to change)

-Represent the dependence of the sinusoids parameters with the energy

It is necesary to install the *[Stingray][1]* Python package.

[1]: https://stingray.readthedocs.io/en/latest/

## binary_cor.py

Corrects the time values due to the binary orbit delay.

## pulse_profile.py

Calculate the pulse profile and perform the fit by sinusoids.

##  XMM_data_selector.py / NuSTAR_data_selector.py

Selects the time of the data fits and return fits files with the corrected time values using binary_cor.py.

XMM_data_selector also selects the XMM-Newton PI channel values and return a fits file with this for a faster use by another programs.

This programs plot the Z<sup>2</sup> and the EF statistics variation with a range of frequencies to search the correct frequency.

##  XMM_pulse_profile.py / NuSTAR_pulse_profile.py

Chooses photons by energy ranges and plot their pulse profiles using pulse_profile.py.

Makes an adjustment of the pulse profiles with sinusoids and return a fits file with their parameters and uncertainties.

##  parameters.py

Represents the variation of the parameters with the energy.
