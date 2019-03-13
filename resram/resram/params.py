import numpy as np
import sys
import boltz_states as bs

## INPUT AND CONSTANTS ##

paramsIn = sys.argv[1]
deltasIn = sys.argv[2]
freqsIn = sys.argv[3]
pumpsIn = sys.argv[4]
		
wg = np.asarray(np.loadtxt(str(freqsIn))) # Ground state normal mode frequencies cm^-1 
we = np.asarray(np.loadtxt(str(freqsIn))) # Excited state normal mode frequencies cm^-1
delta = np.asarray(np.loadtxt(str(deltasIn))) # Dimensionless displacements 
S = (delta**2)/2

with open(paramsIn,'r') as i:
	
	inp = i.readlines()

	j=0
	for l in inp:
		l = l.partition('#')[0]
		l = l.rstrip()
		inp[j] = l
		j+=1
	
	hbar =  5.3088 # plancks constant cm^-1*ps
	T = float(inp[13]) # Temperature K
	kbT = 0.695*T # kbT energy (cm^-1/K)*cm^-1=cm^-1
	cutoff = kbT*0.1 # cutoff for boltzmann dist in wavenumbers
	if T > 10.0:
		beta = 1/kbT # beta cm
		eta = 1/(np.exp(wg/kbT)-1) # array of average thermal occupation numbers for each mode
	elif T < 10.0:
		beta = 1/kbT
		#beta = float("inf")
		eta = np.zeros(len(wg))

	gamma = float(inp[0]) # Homogeneous broadening parameter cm^-1
	theta = float(inp[1]) # Static inhomogenous broadening parameter cm^-1
	E0 = float(inp[2]) # E0 cm^-1

	## Brownian Oscillator parameters ##
	k = float(inp[3]) # kappa parameter
	D =  gamma*(1+0.85*k+0.88*k**2)/(2.355+1.76*k) # D parameter 
	L =  k*D # LAMBDA parameter

	s_reorg = beta*(L/k)**2/2 # reorganization energy cm^-1
	w_reorg = 0.5*np.sum((delta)**2*wg) # internal reorganization energy
	reorg =  w_reorg + s_reorg # Total reorganization energy

	## Time and energy range stuff ##
	ts = float(inp[4]) # Time step (ps)
	ntime = float(inp[5]) #175 # ntime steps
	UB_time = ntime*ts # Upper bound in time range
	t = np.linspace(0,UB_time,ntime) # time range in ps
	EL_reach = float(inp[6]) # How far plus and minus E0 you want 
	EL = np.linspace(E0-EL_reach,E0+EL_reach,1000) # range for spectra cm^-1	
	E0_range = np.linspace(-EL_reach*0.5,EL_reach*0.5,501)# static inhomogeneous convolution range	

	th = np.array(t/hbar) # t/hbar 

	ntime_rot = ntime/np.sqrt(2)
	ts_rot = ts/np.sqrt(2)
	UB_time_rot = ntime_rot*ts_rot
	tp = np.linspace(0,UB_time_rot,ntime_rot) 
	
	tm = np.append(-np.flip(tp[1:],axis=0),tp)
	convEL = np.linspace(E0-EL_reach*0.5,E0+EL_reach*0.5,(max(len(E0_range),len(EL))-min(len(E0_range),len(EL))+1)) # Excitation axis after convolution with inhomogeneous distribution

	M = float(inp[7]) # Transition dipole length angstroms
	n = float(inp[8]) # Refractive index

	rpumps = np.loadtxt(pumpsIn) # Raman pump wavelengths to compute spectra at
	rshift = np.arange(float(inp[9]),float(inp[10]),float(inp[11])) # range and step size of Raman spectrum
	res = float(inp[12]) # Peak width in Raman spectra

	# Determine order from Boltzmann distribution of possible initial states #
	convergence = float(inp[14]) # desired boltzmann coefficient for cutoff
	boltz_toggle = int(inp[15])

	if boltz_toggle == 1:
		boltz_states,boltz_coef,dos_energy = bs.boltz_states()
		if T == 0.0:
			state = 0
		else:
			state = min(range(len(boltz_coef)),key=lambda j:abs(boltz_coef[j]-convergence))

		if state == 0:
			order = 1
		else:
			order = max(max(boltz_states[:state])) + 1
	if boltz_toggle == 0:
		boltz_states,boltz_coef,dos_energy = [0,0,0]
		order = 1

	a = np.arange(order)
	b = a
	Q = np.identity(len(wg),dtype=int)
	wq = np.append(wg,wg)
i.close()


## Prefactors for absorption and Raman cross-sections ##
if order == 1:
	preR = 2.08e-20*(ts**2) #(0.3/pi) puts it in differential cross section 
elif order > 1:
	preR = 2.08e-20*(ts_rot**2)

preA = ((5.744e-3)/n)*ts
preF = preA*n**2

