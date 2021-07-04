def delete_space(A):
	while A[0]==' ':
		if len(A)==1:
			break
		A=A[1:len(A)]
	while A[len(A)-1]==' ':
		if len(A)==1:
			break
		A=A[0:len(A)-1]
	return(A)
def add_space(A):
	word=''
	for j in range(len(A)):
		if A[j] !=' ':word=word+A[j]
		else:word=word+'_'
	return(word)
def lister(A):
	lista=[]
	word=''
	for j in range(len(A)):
		if A[j] ==',':
			lista.append(word)
			word=''
		else:
			word=word+A[j]
		if j==len(A)-1:
			lista.append(word)
			word=''
	return(lista)
def read_files():
	
	try:
		CMD = open('run.sh')
	except:
		CMD = open('run.cmd')
	lista=[]
	run_file=[]

	for linea in CMD:
		if linea[0]==' ':
			continue
		else:
			medidor=0
			cuenta=0
			numero=''
			for k in range(len(linea)):
				if linea[0]==' ' or linea[0]=='#' or linea[0:3]=='rem':
					continue
				if linea[k]!='[' and linea[k]!=']' or medidor==2:
					numero=numero+linea[k]
				if linea[k]=='=':
					lista.append(numero)
					numero=''
					medidor=1
				elif linea[k:k+3]=='rem' or linea[k]=='#' or k==len(linea)-1:
					if medidor==1:
						lista.append(numero[0:len(numero)-1])
						numero=''
						run_file.append(lista)
						medidor=2
				else:
					continue
		
		lista=[]
	
	for j in range(len(run_file)):
		if run_file[j][0][0:8]=='XMM_file':
			XMM_file=lister(delete_space(run_file[j][1]))
		elif run_file[j][0][0:11]=='NuSTAR_file':
			NuSTAR_file=lister(delete_space(run_file[j][1]))
		elif run_file[j][0][0:6]=='source':
			source=delete_space(run_file[j][1])
		elif run_file[j][0][0:5]=='asini':
			asini=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:4]=='Porb':
			Porb=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:3]=='ecc':
			ecc=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:5]=='omega':
			omega=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:2]=='T0':
			T0=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:6]=='period' and run_file[j][0][6]!='_':
			period=0 #float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:4]=='nbin':
			nbin=int(delete_space(run_file[j][1]))
		elif run_file[j][0][0:10]=='nsinusoids':
			nsinusoids=int(delete_space(run_file[j][1]))
		elif run_file[j][0][0:19]=='pulse_frequency_XMM':
			pulse_frequency_XMM=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:22]=='pulse_frequency_NuSTAR':
			pulse_frequency_NuSTAR=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:13]=='period_ranges':
			period_ranges=lister(delete_space(run_file[j][1]))
		elif run_file[j][0][0:11]=='period_bins':
			period_bins=float(delete_space(run_file[j][1]))		
		elif run_file[j][0][0:10]=='XMM_Tstart':
			XMM_time=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:13]=='NuSTAR_Tstart':
			nustar_time=float(delete_space(run_file[j][1]))
		elif run_file[j][0][0:9]=='overwrite':
			overwrite=delete_space(run_file[j][1])
			if overwrite=='N' or overwrite=='n':
				overwrite='n'
			else:
				overwrite='Y'
		elif run_file[j][0][0:9]=='Z_2_check':
			Z_2_check=delete_space(run_file[j][1])
			if Z_2_check=='Y' or Z_2_check=='y':
				Z_2_check='Y'
			else:
				Z_2_check='N'
		elif run_file[j][0][0:13]=='energy_ranges':
			energy_ranges=lister(delete_space(run_file[j][1]))
	period=0
	return(NuSTAR_file,XMM_file,source,asini,Porb,ecc,omega,T0,period,pulse_frequency_NuSTAR,pulse_frequency_XMM,
	energy_ranges,nbin,nsinusoids,overwrite,Z_2_check,period_ranges,period_bins,XMM_time,nustar_time)

def write_files(word,value):
	
	file_to_write=[]
	
	try:
		CMD = open('run.sh')
		key=1
	except:
		CMD = open('run.cmd')
		key=2
	lista=[]
	run_file=[]

	for linea in CMD:
		if linea[0:len(word)]==word:
			pre_linea=linea
			post_linea=''
			for k in range(len(linea)-2):
				if linea[k:k+3]=='rem' or linea[k:k+3]=='REM' or linea[k]=='#':
					post_linea=linea[k:len(linea)]
					pre_linea=linea[0:k]
					break
			for j in range(len(pre_linea)):
				if pre_linea[j]=='=':
					pre_linea=pre_linea[0:j+1]+str(value)+' '
					break
			linea=pre_linea+post_linea
		file_to_write.append([linea])
	if key==1:
		CMD2 = open('run.sh','w')
	if key==2:
		CMD2 = open('run.cmd','w')
	for j in range(len(file_to_write)):
		
		CMD2.write(file_to_write[j][0])
	return()
