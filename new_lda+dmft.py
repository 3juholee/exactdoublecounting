from scipy import *
from h2utils import *
from scipy import linalg
from scipy import optimize
import sys
from excor import ExchangeCorrelation
from scipy import integrate
import time
from scipy import weave
from include_DMFT import *

def DensityLocal(Nc, ist):
	Nc_array = zeros((len(Sol),len(Sol)),dtype = 'float')
	for i,(i1,i2) in enumerate(ist):
		Nc_array[i1,i2] = Nc[i]

	N_local1 = 0.5*(Nc_array[0,0]+Nc_array[0,1]+Nc_array[1,0]+Nc_array[1,1])
	N_local2 = 0.5*(Nc_array[0,0]-Nc_array[0,1]-Nc_array[1,0]+Nc_array[1,1])

	return (N_local1,N_local2)



def Cmp_XC_Potential(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,FX=1.0,FC=1.0):
	"Computes XC potential within LDA"
	code1="""
	        #line 16 "LDA.py"
	        double rho=0.0;
	        for (int i=0; i<an12.extent(0); i++){
	            int i1 = andex(an12(i,0));
	            int i2 = andex(an12(i,1));
	            rho += MM(i1,ix)*LL(i1,ir)*Nc(i)*MM(i2,ix)*LL(i2,ir);
	        }
	        return_val = rho/(2*M_PI*R0*R0*R0/8);
	"""
	code2="""
	    #line 26 "LDA.py"
	    for (int ix=0; ix<mumesh.size(); ix++){
	        for (int ir=0; ir<Ximesh.size(); ir++){
	            vint_xc(ix,ir) = MM(j1,ix)*LL(j1,ir)*MM(j2,ix)*LL(j2,ir)*Vxc(ix,ir)*(Ximesh(ir)*Ximesh(ir)-mumesh(ix)*mumesh(ix));
	        }
	    }
	"""
	def Function1(ir, ix, Nc, in12, index, MM, LL ):
	    rho=0.0
	    for i, (i1,i2) in enumerate(in12):
	        j1 = index[i1]
	        j2 = index[i2]
	        rho += MM[j1,ix]*LL[j1,ir]*Nc[i]*MM[j2,ix]*LL[j2,ir]
	    rs = pow(3/(4*pi*rho),1/3.)
	    return (rs,rho)
	
	def Function2(j1,j2,mumesh,Ximesh,Vxc,MM,LL):
	    vint_xc = zeros((len(mumesh),len(Ximesh)),dtype=float)
	    for ix,eta in enumerate(mumesh):
	        for ir,xi in enumerate(Ximesh):
	            vint_xc[ix,ir] = MM[j1,ix]*LL[j1,ir]*MM[j2,ix]*LL[j2,ir]*Vxc[ix,ir]*(xi**2-eta**2)
	    return vint_xc
	
	exc = ExchangeCorrelation()
	Vxc = zeros((len(mumesh),len(Ximesh)),dtype=float)
	
	an12 = array(in12)
	andex = array(index)
	rho=0.0
	rhox = zeros((len(Ximesh),len(mumesh)),dtype=float)
	for ix in range(len(mumesh)):
	    for ir in range(len(Ximesh)):
	        rho = weave.inline(code1, ['ir','ix','Nc','an12','andex','MM','LL','R0'], type_converters=weave.converters.blitz, compiler = 'gcc')
	        if rho<1e-300:
	            print 'WARNING: rho<0 for ', mumesh[ix], Ximesh[ir], 'rho=', rho
	            rho=1e-300
	        rs = pow(3/(4*pi*rho),1/3.)
	
	        Vxc[ix,ir] = 2*exc.Vx(rs)*FX + 2*exc.Vc(rs)*FC
	        rhox[ir,ix]=rho*(Ximesh[ir]**2-mumesh[ix]**2)*(2*pi*R0**3/8)
	
	rhoz = zeros(len(Ximesh),dtype=float)
	for ir in range(len(Ximesh)):
	    rhoz[ir] = integrate.simps(rhox[ir,:],mumesh)
	print 'Integral(rho)=', integrate.simps(rhoz,Ximesh)
	
	Vxc2 = zeros((len(Sol),len(Sol)))
	for j1 in range(len(Sol)):
		for j2 in range(j1,len(Sol)):
			if Modd[j1]==Modd[j2]:
				vint_xc = zeros((len(mumesh),len(Ximesh)),dtype=float)
				weave.inline(code2, ['vint_xc','j1','j2','mumesh','Ximesh','Vxc','MM','LL','R0'], type_converters=weave.converters.blitz, compiler = 'gcc')
				vint_xc2 = zeros(len(mumesh),dtype=float)
				for ix,eta in enumerate(mumesh):
				    vint_xc2[ix] = integrate.romb(vint_xc[ix,:]*(Ximesh-1.),dhXi)
				Vxc2[j1,j2] = integrate.romb(vint_xc2,mumesh[1]-mumesh[0])
	
	Vxc_ = zeros(len(in12),dtype=float)
	for i, (i1,i2) in enumerate(in12):
	    j1 = index[i1]
	    j2 = index[i2]
	    if (j2>=j1): Vxc_[i] = Vxc2[j1,j2]
	    else: Vxc_[i] = Vxc2[j2,j1]
	
	return Vxc_

def Cmp_XC_Energy(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,FX=1.0,FC=1.0):
    "Computes XC energy within LDA"
    code1="""
            #line 117 "LDA.py"
            double rho=0.0;
            for (int i=0; i<an12.extent(0); i++){
                int i1 = andex(an12(i,0));
                int i2 = andex(an12(i,1));
                rho += MM(i1,ix)*LL(i1,ir)*Nc(i)*MM(i2,ix)*LL(i2,ir);
            }
            return_val = rho/(2*M_PI*R0*R0*R0/8);
    """
    exc = ExchangeCorrelation()
    Exc = zeros((len(mumesh),len(Ximesh)),dtype=float)

    an12 = array(in12)
    andex = array(index)

    for ix,eta in enumerate(mumesh):
        for ir,xi in enumerate(Ximesh):
            rho = weave.inline(code1, ['ir', 'ix', 'Nc','an12','andex','MM','LL','R0'], type_converters=weave.converters.blitz, compiler = 'gcc')
            if rho<0.0:
                print 'WARNING: rho<0 for ', eta, xi, 'rho=', rho
                #print 'Nc=', Nc
                #dsum=0.0
                #for i,(i1,i2) in enumerate(in12):
                #    print i1, i2, MM[i1,ix]*LL[i1,ir]*Nc[i]*MM[i2,ix]*LL[i2,ir], dsum
                #    dsum += MM[i1,ix]*LL[i1,ir]*Nc[i]*MM[i2,ix]*LL[i2,ir]
                rho=1e-300
                
            rs = pow(3/(4*pi*rho),1/3.)
            Exc[ix,ir] = (2*exc.Ex(rs)*rho*FX + 2*exc.Ec(rs)*rho*FC)*(xi**2-eta**2)*(2*pi*R0**3/8)
            
    Ene = zeros(len(mumesh),dtype=float)
    for ix,eta in enumerate(mumesh):
        Ene[ix] = integrate.romb(Exc[ix,:]*(Ximesh-1.),dhXi)

    return integrate.romb(Ene,mumesh[1]-mumesh[0])



def Cmp_local_XC_Potential(Nc,MM,LL,mumesh,Ximesh,Project_i,Project_f,R0):
	"Computes XC potential within LDA"
	code3="""
		#line 126 "h2_LDA.py"
		double rho=0.0;
		for (int j=0; j<Project_i.extent(0); j++){
		    int i=Project_i(j,0);
		    int i1=Project_i(j,1);
		    int i2 = Project_i(j,2);
		    rho += MM(i1,ix)*LL(i1,ir)*Nc(i)*MM(i2,ix)*LL(i2,ir)*Project_f(j);
		}
	return_val = rho/(2*M_PI*R0*R0*R0/8);
	"""
	code4="""
		#line 137 "h2_LDA.py"
		for (int ix=0; ix<mumesh.size(); ix++){
		    for (int ir=0; ir<Ximesh.size(); ir++){
		        vint_x(ix,ir) = MM(j1,ix)*LL(j1,ir)*MM(j2,ix)*LL(j2,ir)*Vx(ix,ir)*(Ximesh(ir)*Ximesh(ir)-mumesh(ix)*mumesh(ix));
		        vint_c(ix,ir) = MM(j1,ix)*LL(j1,ir)*MM(j2,ix)*LL(j2,ir)*Vc(ix,ir)*(Ximesh(ir)*Ximesh(ir)-mumesh(ix)*mumesh(ix));
		    }
		}
	"""
	exc = ExchangeCorrelation()
	Vx = zeros((len(mumesh),len(Ximesh)),dtype=float)
	Vc = zeros((len(mumesh),len(Ximesh)),dtype=float)
	
	N_local1,N_local2 = DensityLocal(Nc,in12)

	rho=0.0
	rhox = zeros((len(mumesh),len(Ximesh)),dtype=float)
	for ix in range(len(mumesh)):
		for ir in range(len(Ximesh)):
			#rho = N_LL*(1/sqrt(2)(phi(1)+phi(2)))**2
			rho = N_local1*0.5*(MM[0,ix]*LL[0,ir]+MM[1,ix]*LL[1,ir])**2
			rho *= 1/(2*pi*R0*R0*R0/8) 
			#(R/2)^3 comes from h2_iter.py
			#2*pi comes from phi = 1/sqrt(2pi)LL(xi)MM(eta)


			if rho<0.0:
			    print 'WARNING: rho<0 for ', ix, ir, 'rho=', rho
			    rho=1e-300
			
			rs = pow(3/(4*pi*rho),1/3.)
			
			Vx[ix,ir] = 2*exc.Vx(rs)
			Vc[ix,ir] = 2*exc.Vc(rs)
			rhox[ix,ir]=rho*(Ximesh[ir]**2-mumesh[ix]**2)*(2*pi*R0**3/8)

		
	rhoz = zeros(len(Ximesh),dtype=float)
	for ir in range(len(Ximesh)):
		rhoz[ir] = integrate.simps(rhox[:,ir],mumesh)
	print 'Integral(rho_local)=', integrate.simps(rhoz,Ximesh)
	
	Vx_local=0.0
	Vc_local=0.0

	vint_x = rhox/N_local1*Vx
	vint_c = rhox/N_local1*Vc
	vint_x2 = zeros(len(mumesh),dtype=float)
	vint_c2 = zeros(len(mumesh),dtype=float)
	for ix,eta in enumerate(mumesh):
		vint_x2[ix] = integrate.romb(vint_x[ix,:]*(Ximesh-1.),dhXi)
		vint_c2[ix] = integrate.romb(vint_c[ix,:]*(Ximesh-1.),dhXi)
	
	Vx_local = integrate.romb(vint_x2,mumesh[1]-mumesh[0])
	Vc_local = integrate.romb(vint_c2,mumesh[1]-mumesh[0])
	
	return (Vx_local,Vc_local)


def Cmp_local_XC_Energy(Nc,MM,LL,mumesh,Ximesh,Project_i,Project_f,R0):
	"Computes XC energy within LDA"
	code3="""
	        #line 200 "h2_LDA.py"
	        double rho=0.0;
	        for (int j=0; j<Project_i.extent(0); j++){
	            int i=Project_i(j,0);
	            int i1=Project_i(j,1); 
							int i2 = Project_i(j,2);
	            rho += MM(i1,ix)*LL(i1,ir)*Nc(i)*MM(i2,ix)*LL(i2,ir)*Project_f(j);
	        }
	        return_val = rho/(2*M_PI*R0*R0*R0/8);
	"""
	code4="""
	    #line 211 "h2_LDA.py"
	    for (int ix=0; ix<mumesh.size(); ix++){
	        for (int ir=0; ir<Ximesh.size(); ir++){
	            E_x(ix,ir) = MM(j1,ix)*LL(j1,ir)*MM(j2,ix)*LL(j2,ir)*Ex(ix,ir)*(Ximesh(ir)*Ximesh(ir)-mumesh(ix)*mumesh(ix));
	            E_c(ix,ir) = MM(j1,ix)*LL(j1,ir)*MM(j2,ix)*LL(j2,ir)*Ec(ix,ir)*(Ximesh(ir)*Ximesh(ir)-mumesh(ix)*mumesh(ix));
	        }
	    }
	"""
	N_local1,N_local2 = DensityLocal(Nc,in12)
	exc = ExchangeCorrelation()
	Ex = zeros((len(mumesh),len(Ximesh)),dtype=float)
	Ec = zeros((len(mumesh),len(Ximesh)),dtype=float)
	rhox = zeros((len(mumesh),len(Ximesh)),dtype=float)
	for ix,eta in enumerate(mumesh):
		for ir,xi in enumerate(Ximesh):
			#rho = weave.inline(code3, ['ir','ix','Nc','MM','LL','R0','Project_i','Project_f'], type_converters=weave.converters.blitz, compiler = 'gcc')
			
			rho = N_local1*0.5*(MM[0,ix]*LL[0,ir]+MM[1,ix]*LL[1,ir])**2
			rho *= 1/(2*pi*R0*R0*R0/8) 
			#(R/2)^3 comes from h2_iter.py
			#2*pi comes from phi = 1/sqrt(2pi)LL(xi)MM(eta)


			if rho<0.0:
			    print 'WARNING: rho<0 for ', ix, ir, 'rho=', rho
			    rho=1e-300
			
			rs = pow(3/(4*pi*rho),1/3.)
			
			Ex[ix,ir] = 2*exc.Ex(rs)
			Ec[ix,ir] = 2*exc.Ec(rs)
			
			rhox[ix,ir]=rho*(Ximesh[ir]**2-mumesh[ix]**2)*(2*pi*R0**3/8)

		
	rhoz = zeros(len(Ximesh),dtype=float)
	for ir in range(len(Ximesh)):
		rhoz[ir] = integrate.simps(rhox[:,ir],mumesh)
	#print 'Integral(rho)=', integrate.simps(rhoz,Ximesh)
	
	Eint_x = rhox*Ex
	Eint_c = rhox*Ec
	Eint_x2 = zeros(len(mumesh),dtype=float)
	Eint_c2 = zeros(len(mumesh),dtype=float)
	for ix,eta in enumerate(mumesh):
		Eint_x2[ix] = integrate.romb(Eint_x[ix,:]*(Ximesh-1.),dhXi)
		Eint_c2[ix] = integrate.romb(Eint_c[ix,:]*(Ximesh-1.),dhXi)
	
	Ex_local = integrate.romb(Eint_x2,mumesh[1]-mumesh[0])
	Ec_local = integrate.romb(Eint_c2,mumesh[1]-mumesh[0])
	
	return (Ex_local,Ec_local)


if __name__ == '__main__':
	
	if len(sys.argv)<=1:
	    print 'Give input integer iR!'
	    sys.exit(1)
	
	DMFT=True
	RUUN = True
	
	    
	iR=int(sys.argv[1])
	print 'iR=', iR
	
	# DMFT parameters
	beta = 100         # Temperature
	Nomega = 50*beta   # Number of all Matsubara points
	nom = 150          # Number of frequency points to sample
	ntail = 100        # Number of frequency points in the tail
	gamma = 0.0        # broadening
	Noccup=2.0         # occupancy for H2
	mu_a, mu_b = -2., 2. # limits where to look for the chemical potential
	
	
	mixr=0.5     # mixing of the charge density
	Nr_xc=2**8+1  # number of points for integration of V_xc in r-direction
	Nc_xc=2**6+1  # number of points for integration of V_xc in eta-direction
	rmin=1e-7     # mesh-offset
	rmax=20       # maximum r-value in atomic units
	Nitt=100      # maximum number of iterations
	Nitt=50     # maximum number of iterations
	justHF=False  # For testing, we can switch to Hartree-Fock
	
	# Creates necessary real space meshes
	Ximesh = 1 + logspace(log10(rmin),log10(rmax),Nr_xc)
	mumesh = linspace(-1+1e-7,1-1e-7,Nc_xc)
	cc=log(Ximesh-1)
	dhXi = (cc[-1]-cc[0])/(len(cc)-1)  # For romberg integration
	
	# Creates necessary frequency meshes
	omega_large = linspace(pi/beta,pi*(2*Nomega-1)/beta,Nomega)
	print 'Matsubara mesh start end:', omega_large[0], omega_large[-1]
	sigdata = zeros((3,len(omega_large)),dtype=float)
	sigdata[0,:] = omega_large
	
	ctqmc = IMP(nom,beta,RUUN)
	
	# Reads H2+ solutions
	(Sol,Ck,ak) = ReadSolution(iR)
	execfile('H2+/H2_Coulomb.'+str(iR))
	
	
	# Compues local U, i.e., F0
	U_local=0.0
	for i1 in range(2):
	    for i2 in range(2):
	        for i3 in range(2):
	            for i4 in range(2):
	                U_local += Uc[i1,i2,i3,i4]*0.25
	print "U_local = ",U_local
	
	if os.path.isfile(ctqmc.dir+'/'+ctqmc.Signame):
	    # use the old file, if exists
	    sold = loadtxt(ctqmc.dir+'/'+ctqmc.Signame).transpose()
	    for i in range(1,len(sold)):
	        Sf = interpolate.UnivariateSpline(sold[0,:],sold[i,:],s=0)
	        sigdata[i,:] = Sf(omega_large)
	    #s_oo=[float(loadtxt(ctqmc.dir+'/'+'s_oo'))]
	else:
	    sigdata[2,:] = -gamma
	
	if os.path.isfile(ctqmc.dir+'/'+'s_oo'):
	    (s2,Edc) = loadtxt(ctqmc.dir+'/'+'s_oo')
	    s_oo=[s2]
	else:
	    s_oo = [0.5*U_local]
	    Edc = 0.5*U_local
	    
	Sigs,oms_equal,ind_om = create_log_mesh(sigdata, nom, ntail)
	
	index=[]
	for i,(R,Ene,m,p,A) in enumerate(Sol):
	    for deg in range(2*m+1): index.append(i)
	R0 = R
	
	print 'All states size=', len(index), 'index=', index
		
	# Which states are even and which are odd?
	Modd=zeros(len(Sol),dtype=int)
	small=1e-5
	MM = zeros((len(Sol),len(mumesh)),dtype=float)
	LL = zeros((len(Sol),len(Ximesh)),dtype=float)
	for i,(R,Ene,m,p,A) in enumerate(Sol):
	    for ir,r in enumerate(Ximesh): LL[i,ir] = Lambda(r,m,p,Ck[i])
	    for ix,x in enumerate(mumesh): MM[i,ix] = Mm(x,m,p,ak[i])
	    ## Checking which solution is even and which is odd (gerade/ungerade)
	    MMx = [MM[i,ix]*x for ix,x in enumerate(mumesh)]
	    IMM = integrate.romb(MM[i,:],mumesh[1]-mumesh[0])
	    IMMx = integrate.romb(MMx,mumesh[1]-mumesh[0])
	    if abs(IMM*IMMx)>small:
	        print 'ERROR: function M has to have a well defined parity! It does not look like it has for m=', m, 'Ene=', Ene
	        print '<M>=', IMM, '<M*x>=', IMMx
	    if abs(IMM)<small: Modd[i]=1
	print 'Which wave functions are odd =', Modd
	

	# Making index for two states
	ist=[]
	in12=[]
	in12_1={}
	ll=0
	for i,i1 in enumerate(index):
	    for j,i2 in enumerate(index):
	        ist.append( (i1,i2) )
	        in12.append( (i,j) )
	        in12_1[(i,j)]=ll
	        ll+=1
	
	# From even-odd basis to local basis
	T=matrix([[1/sqrt(2.),1/sqrt(2.)],[1/sqrt(2.),-1/sqrt(2.)]])
	
	H0 = zeros((len(index),len(index)),dtype=float)
	V2 = zeros((len(index),len(index)),dtype=float)
	Nt2 = zeros((len(index),len(index)),dtype=float)
	
	for i,ii in enumerate(index):
	    H0[i,i]=Sol[ii][1]-2./R # consider the electron's energy only
	
	# Hartree and Fock potential
	UH = zeros((len(ist),len(ist)), dtype=float)
	UF = zeros((len(ist),len(ist)), dtype=float)
	for i, (i1,i2) in enumerate(ist):
	    for j, (i3,i4) in enumerate(ist):
	        UH[i,j] = Uc[i1,i3,i4,i2]
	        UF[i,j] = Uc[i1,i3,i2,i4]
	
	# Projector
	Project_i=array([[in12_1[(0,0)],0,0],[in12_1[(1,1)],1,1]])
	Project_f=array([0.5,0.5])
	
	# density matrix
	Nc = zeros(len(ist),dtype=float)
	Nc0 = zeros(len(ist),dtype=float)
	
	if os.path.exists(ctqmc.dir+'/'+'Nc.dat'):
	    Nc = loadtxt(ctqmc.dir+'/'+'Nc.dat')
	else:
	    # more appropriate for localized regime
	    Nc[in12_1[(0,0)]]=1.0
	    Nc[in12_1[(1,1)]]=1.0
	
	p_Ener=0
	
	
	
	
	if (True):
	    for i,(i1,i2) in enumerate(in12): Nt2[i1,i2] = Nc[i]
	    
	    print 'Nt2[In H2+ basis]='
	    for i in range(len(Nt2)):
	        for j in range(len(Nt2)):
	            print "%10.5f " % Nt2[i,j].real,
	        print
	
	
	oms = Sigs[0]
	
	FullExchange = True
	if FullExchange:
	    FX=0.0
	    FC=1.0
	else:
	    FX=1.0
	    FC=1.0
	
	
	muf = open('imp/mu.dat','w')
	Nimpf = open('imp/Nimp.dat','w')
	EDMFTf = open('imp/Phi_DMFT.dat','w')
	Eimpf = open('imp/Eimp.dat','w')
	s_oof = open('imp/s_oo.dat','w')
	Etotf = open('imp/Etot.dat','w')
	
	print
	print "start Hartree Fock"
	print
	for itt in range(Nitt):
	
		VHa = dot(UH, Nc)     # Hartree
		VFx = -0.5*dot(UF,Nc) # Fock
		
		Vxc = Cmp_XC_Potential(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,FX,FC) + VFx*(1.-FX)
		
		(Vx_local, Vc_local) = Cmp_local_XC_Potential(Nc,MM,LL,mumesh,Ximesh,Project_i,Project_f,R0)
		(epsx_local, epsc_local) = Cmp_local_XC_Energy(Nc,MM,LL,mumesh,Ximesh,Project_i,Project_f,R0)
		
		for i, (i1,i2) in enumerate(in12): V2[i1,i2]=VHa[i]+Vxc[i]
		
		(Ener,Evec) = linalg.eigh(H0+V2)
		mu=0.5*(Ener[0]+Ener[1])
		
		
		Proj = Evec[0:2,:] # Projector
		Proj = matrix(Proj)
		TProj = T*Proj
		
		# computes Hartree-Fock energy levels
		eks=Cmp_eks(Sigs,s_oo,Edc,Ener,TProj)
		# computes the current chemical potential
		mu = Cmp_mu(Noccup/2.,mu_a,mu_b,eks,oms,oms_equal,beta)
		# computes Gloc & Delta
		(Eimp,Glocal,Delta) = Cmp_Gloc_Delta(Sigs,s_oo,mu,Edc,Ener,TProj)
		 
		
		(Sigs, s_oo, n_imp, TrSigmaG) = ctqmc.Run(U_local,-Eimp,omega_large,ind_om,Glocal,Delta,s_oo)
		ctqmc.copyfiles(itt)
		
		print '2: s_oo=', s_oo, 'n_imp=', n_imp
		
		
		(eks,Ars) = Cmp_eVks(Sigs,s_oo,Edc,Ener,TProj)
		# recomputes the current chemical potential for charge neutrality
		#mu = Cmp_mu(Noccup/2.,mu_a,mu_b,eks,oms,oms_equal,beta)
		#print 'mu[1]=', mu
		
		# In hartree-fock basis
		Nt = 2.*Cmp_DensityMatrix(eks,Ars,mu,oms,oms_equal,beta)  # 2 due to spin
		Nt2 = matrix(Evec) * Nt * matrix(Evec).T                  # From HF to H2+ basis
		for i, (i1,i2) in enumerate(in12):  Nc0[i] = real(Nt2[i1,i2])
		
		#print 'imag(Nc)=', sum(sum(abs(imag(Nt2))))
		 
		# mixing of density
		Nc[:] = mixr*Nc0[:] + (1.-mixr)*Nc[:]
		
		N_local = 0.5*(Nc[in12_1[(0,0)]]+Nc[in12_1[(0,1)]]+Nc[in12_1[(1,0)]]+Nc[in12_1[(1,1)]])
		
		Nf = n_imp
		#Nf = N_local
		Edc   = U_local*Nf/2. + Vc_local
		#PhiDC_local = U_local*Nf*Nf/4. *2 + epsc_local*Nf *2
		PhiDC_local = U_local*Nf*Nf/4. *2 + epsc_local *2
		PhiDMFT = 2*TrSigmaG
		
		
		savetxt(ctqmc.dir+'/'+'s_oo',[s_oo,[Edc]*len(s_oo)])
		
		# Kinetic energy part
		Etot0 = sum([H0[i,i]*real(Nt2[i,i]) for i in range(len(H0))])
		# Hartree Energy
		EHa = 0.5*sum(Nc*VHa)
		# Full Fock Energy
		EFx = 0.5*sum(Nc*VFx)
		# Exchange-correlation energy
		Exc = Cmp_XC_Energy(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,FX,FC)
		
		# Total energy
		Etot = Etot0 + EHa + Exc + (1.-FX)*EFx + PhiDMFT - PhiDC_local + 2./R
		
		Ec = Cmp_XC_Energy(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,0.0,1.0)
		Ex = Cmp_XC_Energy(Nc,MM,LL,mumesh,Ximesh,in12,index,R0,1.0,0.0)
		
		#print 'E0=', Etot0, 'EH=', EHa, 'EFexact=', EFx, 'Ex=', Ex, 'Ec=', Ec, 'Exc=', Exc, 'EDMFT=', PhiDMFT, 'PhiDC=', PhiDC_local
		#print 'Etot=', Etot, 'Ediff=', abs(Etot-p_Ener)
		#
		#print 'N_local=', N_local, 'N_imp=', n_imp
		#print 'Edc,PhiDC_local=', Edc, PhiDC_local
		#print 'Ec=', Ec, 'Ec_local=', epsc_local*n_imp*2, 'Ex=', Ex, 'Ex_exact=', EFx, 'Ex_DMFT=', -U_local*n_imp**2/2
		
			
		print 
		print 'mu=', mu
		print 'Eimp = ',	Eimp
		print 'Nimp = ',	Nf
		print 'N_local-N_imp=', N_local-n_imp, 'N_local=', N_local, 'N_imp=', n_imp
		print 'T[n] = ', Etot0
		print 'EHF[n] = ', EHa+EFx*(1.-FX)
		print 'Ec[n] = ', Exc
		print 'PhiDMFT = ', PhiDMFT
		print 'PhiDC_local = ', PhiDC_local
		print "Etot = ",	Etot
		print R, Etot, abs(Etot-p_Ener), 'itt = ', itt
		savetxt(ctqmc.dir+'/'+'s_oo',[s_oo,[Edc]*len(s_oo)])

		print >> muf, itt, mu
		print >> Eimpf, itt, Eimp
		print >> s_oof, "s_oo["+str(itt)+"]=", s_oo[0]
		print >> Nimpf, itt, Nf
		print >> EDMFTf, itt, PhiDMFT
		print >> Etotf, itt, Etot
		 
		savetxt(ctqmc.dir+'/'+'Nc.dat',Nc)
		if abs(Etot-p_Ener)<5e-6: break
		p_Ener = Etot
