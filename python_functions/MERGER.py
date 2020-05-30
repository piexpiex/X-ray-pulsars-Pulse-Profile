from PyPDF2 import PdfFileMerger
from read_files import *
READ=read_files()
source=add_space(READ[2])


print('Merging the pdfs with the resulting plots')

pdfs = [source+"_XMM_pulse_profile.pdf",source+"_XMM_pulse_profile_harmonics.pdf","NuSTAR_pulse_profile.pdf","NuSTAR_pulse_profile_harmonics.pdf","amplitudes_energy_variation.pdf"]
pdf_name ='figures_folder/'+ source+'_summary.pdf'
merger = PdfFileMerger()

for pdf in pdfs:
	try:
		merger.append(open('figures_folder/'+pdf, 'rb'))
	except:
		continue

with open(pdf_name, 'wb') as name:
    merger.write(name)
