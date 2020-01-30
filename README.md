# X-ray-pulsars-Pulse-Profile

Programs used to obtain pulse profiles of X-ray pulsar with NuSTAR and XMM-Newton data.

These programs include:

-Select data from Nustar or XMM-Newton fits files

-Search frequency of the pulsar

-Represent the pulse profiles in differents energy ranges

-Adjusment of pulse profiles by five sinusoids

#binary_cor.py

Corrects the time values due to the binary orbit delay.

It is recommended to have installed PyAstronomy for faster run of the program.

#XMM_data_selector.py / NuSTAR_data_selector.py

Selects the time of the data fits and return fits files with the corrected time values

XMM_data_selector also selects the XMM Newton PI channel values and return a fits file with this for a faster use by another programs.
