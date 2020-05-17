from PyPDF2 import PdfFileMerger
from read_files import *
READ=read_files()
source=add_space(READ[2])
print(READ[2])
print(source)

pdfs = [source+"_XMM_pulse_profile.pdf",source+"_XMM_pulse_profile_harmonics.pdf","NuSTAR_pulse_profile.pdf","NuSTAR_pulse_profile_harmonics.pdf","final.pdf"]
nombre_archivo_salida = source+'.pdf'
fusionador = PdfFileMerger()

for pdf in pdfs:
    fusionador.append(open(pdf, 'rb'))

with open(nombre_archivo_salida, 'wb') as salida:
    fusionador.write(salida)
