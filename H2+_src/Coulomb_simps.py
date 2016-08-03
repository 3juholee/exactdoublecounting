#from original h2_Coulomb.py, the romb integration of cmp_Coulomb
#is changed to simps integration.



import sys
from scipy import *
from scipy import special, integrate
from pylab import *
import poly as polyx
from time import time
from scipy.misc import derivative
from scipy import weave


def simps(Y,Xx):
	if len(Xx)!=len(Y): 
		print "the # of samplings of X and Y are different"
		return
	res = 0.0
	
	simps = """
		#line 16 "simps_weave.py"
		int step = 2;
		int N = Xx.size();
		double x1;
		double x2;
		double res = 0.0;
		
		if ( N % 2 != 0){
			for (int i=0; i < Xx.size()-2; i+=step)
			{
				x1 = Xx(i+1)-Xx(i);
				x2 = Xx(i+2)-Xx(i+1);
				res += (x1+x2)/6.*( Y(i)*(2.-x2/x1) + Y(i+1)*((x1+x2)*(x1+x2)/(x1*x2))+Y(i+2)*(2.-x1/x2));
			}
		}
		
		else{
			for (int i=0; i < N-3; i+=step)
			{
				x1 = Xx(i+1)-Xx(i);
				x2 = Xx(i+2)-Xx(i+1);
				res += (x1+x2)/6.*( Y(i)*(2.-x2/x1) + Y(i+1)*((x1+x2)*(x1+x2)/(x1*x2))+Y(i+2)*(2.-x1/x2));
			}
			res += 0.5*(Xx(N-1)-Xx(N-2))*(Y(N-1)+Y(N-2));
		}
		return_val = res;
	"""

	res = weave.inline(simps,['Y','Xx'],  type_converters= weave.converters.blitz, compiler = 'gcc')
	return res


def prolate((x,y,z),R):
    "Converts from cartesian coordinates to prolate"
    r2 = (x**2+y**2+z**2)/(R/2)**2
    z2 = (z/(R/2))**2
    sq=sqrt((1+r2)**2-4*z2)
    
    xi=sqrt(0.5*(1+r2+sq))
    et=sqrt(0.5*(1+r2-sq))*sign(z)
    ph=arctan2(y,x)
    return (xi,et,ph)

def cartesian((xi,eta,ph),R):
    "Converts from prolate to cartesian coordinates"
    rho=0.5*R*sqrt(xi**2-1)*sqrt(1-eta**2)
    return (rho*cos(ph),rho*sin(ph),0.5*R*xi*eta)

def Mm(x,m,q,aks):
    "Function M(eta), part of the H2+ solution"
    m1 = abs(m)
    maxk=len(aks)
    #Plm = polyx.Plmc(x,m,maxk)
    #mm = [aks[il]*Plm[il] for il in range(len(aks))]
    Plm = special.lpmn(m1, maxk-1+m1, x)[0][m1]
    mm = [aks[il]*Plm[il+m1] for il in range(len(aks))]
    return sum(mm)*exp(-q*x)

def Lambda(lam,m,p,Cks):
    "Function Lambda(lambda),part of the H2+ solution"
    m1 = abs(m)
    alpha = m1+0.0
    max_n = m1+len(Cks)-1
    x = 2*p*(lam-1.)
    Lag = polyx.Laguerre_dirty(x,alpha,max_n)
    #csum = sum([Cks[k]*Lag[k+m] for k in range(len(Cks))])  ### original
    csum = sum([Cks[k]*Lag[k] for k in range(len(Cks))])  ### HERE Lag[k]->Lag[m+k]
    return sqrt((lam**2-1)**m1)*exp(-x/2.)*csum

def ReadSolution(iR):
    fs = open('H2+/H2_solution.'+str(iR), 'r')
    fs.readline()
    Sol=[]
    Ck=[]
    ak=[]
    for il,line in enumerate(fs):
        exec(line)
        Sol.append( sol )
        Ck.append( Cks )
        ak.append( aks )
    return (Sol,Ck,ak)

def cmp_MPM(Ny, maxm,maxl, ist, ak, Modd):
    Xx = linspace(-1,1,Ny)
    
    MPM =  zeros((len(ist),maxl+1),dtype=float)
    MPxM = zeros((len(ist),maxl+1),dtype=float) 
    for i,(i1,i2) in enumerate(ist):
        (R1,Ene1,m1,p1,A1) = Sol[i1]
        (R2,Ene2,m2,p2,A2) = Sol[i2]
        
        print 'i1=', i1, 'i2=', i2, 'm1=', m1, 'm2=', m2 #, MPM,MPxM
        m = abs(m1-m2)
        
        Plm_all = array([special.lpmn(maxm,maxl,x)[0] for x in Xx])
        MM = array([Mm(x,m1,p1,ak[i1]) * Mm(x,m2,p2,ak[i2]) for x in Xx])
        Plm = transpose(Plm_all[:,m])
        for il in range(m,maxl+1):
            odd = (Modd[i1]+Modd[i2]+il+m)%2
            if not odd:  # If the function is odd, we do not calculate!
                MMP = MM*Plm[il,:]
                MMPx = MMP * Xx**2
                MPM [i,il] = integrate.romb(MMP, Xx[1]-Xx[0])
                MPxM[i,il] = integrate.romb(MMPx,Xx[1]-Xx[0])
    
    return (MPM,MPxM)

def cmp_LPQL(Nx,xmin,xmax, maxm,maxl, ist, Sol, Ck, Modd):
	
	Xx=1-xmin+logspace(log10(xmin),log10(xmax),Nx)
	# For rombohedral method we need the change of variable to equidistant mesh
	# We can change the variable Int[f(x),{x,1,Inf}] = Int[ exp(t) f(1-xmin+exp(t)), {t, log(xmin),inf}]
	# Here t is distributed in linear mesh with spacing dxi, and exp(t)=x-1+xmin, because we
	# constructued logarithmic mesh by x=1-xmin+exp(t)
	Expu = Xx-1.+xmin 
	dxi = (log(Expu[-1])-log(Expu[0]))/(len(Expu)-1.)
	
	Qlm_all = array([special.lqmn(maxm,maxl,x)[0] for x in Xx])
	Plm_all = array([special.lpmn(maxm,maxl,x)[0] for x in Xx])
	Qlm_all[0]=0.0  # This diverges, hence we set it to zero
	
	LPL = zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LQL = zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LPxL= zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LQxL= zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LLP = zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LLQ = zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LLPx= zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	LLQx= zeros((len(ist),maxl+1,len(Xx)),dtype=float)
	
	for i,(i1,i2) in enumerate(ist):
		(R1,Ene1,m1,p1,A1) = Sol[i1]
		(R2,Ene2,m2,p2,A2) = Sol[i2]
		print 'i1=', i1, 'i2=', i2, 'm1=', m1, 'm2=', m2
		
		m = abs(m1-m2)
		
		Qlm = transpose(Qlm_all[:,m])
		Plm = transpose(Plm_all[:,m])
		LL   = array([Lambda(x,m1,p1,Ck[i1]) * Lambda(x,m2,p2,Ck[i2]) for ix,x in enumerate(Xx)])
		
		for il in range(m,maxl+1):
			odd = (Modd[i1]+Modd[i2]+il+m)%2
			
			if odd: continue
			

			LLQ [i,il,:] = LL * Qlm[il,:]
			LLQx[i,il,:] = LLQ[i,il,:] * Xx**2
			LLP [i,il,:] = LL * Plm[il,:]
			LLPx[i,il,:] = LLP[i,il,:] * Xx**2


			for ix in range(len(Xx)):
				LPL [i,il,ix] = simps(LLP [i,il,:ix+1],Xx[:ix+1])
				LPxL[i,il,ix] = simps(LLPx[i,il,:ix+1],Xx[:ix+1])
				LQL [i,il,ix] = simps(LLQ [i,il,ix:],Xx[ix:])
				LQxL[i,il,ix] = simps(LLQx[i,il,ix:],Xx[ix:])
	    
	
	return (LPL, LQL, LPxL, LQxL, LLP, LLQ, LLPx, LLQx)


def cmp_Coulomb(Nx,xmin,xmax, maxm,maxl, ist, Sol, Modd, MPM, MPxM, LPL, LQL, LPxL, LQxL, LLP, LLQ, LLPx, LLQx):
    
    Xx=1-xmin+logspace(log10(xmin),log10(xmax),Nx)
    # For rombohedral method we need the change of variable to equidistant mesh
    # We can change the variable Int[f(x),{x,1,Inf}] = Int[ exp(t) f(1-xmin+exp(t)), {t, log(xmin),inf}]
    # Here t is distributed in linear mesh with spacing dxi, and exp(t)=x-1+xmin, because we
    # constructued logarithmic mesh by x=1-xmin+exp(t)
    Expu = Xx-1.+xmin 
    dxi = (log(Expu[-1])-log(Expu[0]))/(len(Expu)-1.)
    
    Uc = zeros((len(Sol),len(Sol),len(Sol),len(Sol)),dtype=float)
    for i,(i1,i2) in enumerate(ist):
        for j,(i3,i4) in enumerate(ist):
            (R1,Ene1,m1,p1,A1) = Sol[i1]
            (R2,Ene2,m2,p2,A2) = Sol[i2]
            (R3,Ene3,m3,p3,A3) = Sol[i3]
            (R4,Ene4,m4,p4,A4) = Sol[i4]
            if (m1-m2!=m4-m3): continue
            m = abs(m1-m2)
            
            uc=0
            for il in range(m,maxl+1):
                odd1 = (Modd[i1]+Modd[i2]+il+m)%2
                odd2 = (Modd[i3]+Modd[i4]+il+m)%2
                if odd1 or odd2: continue
                
                r1 = MPM [i,il]*MPM [j,il] * (simps(LLQx[i][il] * LPxL[j][il],Xx) + simps(LLPx[i][il] * LQxL[j][il],Xx))
                r2 = MPxM[i,il]*MPxM[j,il] * (simps(LLQ [i][il] * LPL [j][il],Xx) + simps(LLP [i][il] * LQL [j][il],Xx))
                r3 = MPM [i,il]*MPxM[j,il] * (simps(LLQx[i][il] * LPL [j][il],Xx) + simps(LLPx[i][il] * LQL [j][il],Xx))
                r4 = MPxM[i,il]*MPM [j,il] * (simps(LLQ [i][il] * LPxL[j][il],Xx) + simps(LLP [i][il] * LQxL[j][il],Xx))

                res = (-1)**m * (2*il+1)*(math.factorial(il-m)/(math.factorial(il+m)+0.0))**2*(r1+r2-r3-r4)
                uc += res
            #print  'i1=', i1, 'i2=', i2, 'i3=', i3, 'i4=', i4, Uc[i1,i2,i3,i4]  
            Uc[i1,i3,i4,i2] = uc
            
    (R,Ene0,m0,p0,A0) = Sol[0]
    Uc *= 4./R

    return Uc





#################start ################################

execfile('./input.py')


xmin=1e-7
xmax=20
Nx=2**9+1
Ny=2**8+1

maxl=5
maxm=3


execfile('./which.py')
for iR in which:
	print 'iR=', iR, 'maxl=', maxl, 'maxm=', maxm, 'Nx=', Nx, 'Ny=', Ny
	
	(Sol,Ck,ak) = ReadSolution(iR)
	
	pmin=1e100
	for il in range(len(Sol)):
	    (R,Ene,m,p,A) = Sol[il]
	    pmin = min(p,pmin)
	xmax = max(xmax,30./pmin)
	print 'xmax=', xmax
	
	# Checking which solution is even and which is odd (gerade/ungerade)
	Modd=zeros(len(Sol),dtype=int)
	small=1e-8
	for i,sol in enumerate(Sol):
	    (R,Ene,m,p,A) = sol
	    Xx = linspace(-1,1,Ny)
	    MM = [Mm(x,m,p,ak[i]) for x in Xx]
	    MMx = [MM[ix]*x for ix,x in enumerate(Xx)]
	    IMM = integrate.romb(MM,Xx[1]-Xx[0])
	    IMMx = integrate.romb(MMx,Xx[1]-Xx[0])
	    if abs(IMM*IMMx)>small:
	        print 'ERROR: function M has to have a well defined parity! It does not look like it has for m=', m, 'Ene=', Ene
	    if abs(IMM)<small: Modd[i]=1
	print 'Which wave function is odd =', Modd
	
	# Making index for two states
	ist=[]
	for i1 in range(len(Sol)):
	    for i2 in range(len(Sol)):
	        ist.append( (i1,i2) )
	
	# Computing integrals for angular part <M|P|M>
	(MPM,MPxM) = cmp_MPM(Ny, maxm,maxl, ist, ak, Modd)
	# Computing integral for the radial part <Lambda|P|Lambda>
	(LPL, LQL, LPxL, LQxL, LLP, LLQ, LLPx, LLQx) = cmp_LPQL(Nx,xmin,xmax, maxm,maxl, ist, Sol, Ck, Modd)
	# Combining angular and radial part
	Uc = cmp_Coulomb(Nx,xmin,xmax, maxm,maxl, ist, Sol, Modd, MPM, MPxM, LPL, LQL, LPxL, LQxL, LLP, LLQ, LLPx, LLQx)
	
	
	exact=2
	# Below is just printing
	print 'U which needs exact treatment'
	for i,(i1,i2) in enumerate(ist):
	    (R1,Ene1,m1,p1,A1) = Sol[i1]
	    (R2,Ene2,m2,p2,A2) = Sol[i2]
	    for j,(i3,i4) in enumerate(ist):
	        (R3,Ene3,m3,p3,A3) = Sol[i3]
	        (R4,Ene4,m4,p4,A4) = Sol[i4]
	        if i1<exact and i2<exact and i3<exact and i4<exact:
	            if abs(Uc[i1,i2,i3,i4])>1e-7:
	                print '(', i1, ',', i2, ',', i3, ',', i4, ') ->', Uc[i1,i2,i3,i4]
	
	fo = open('H2+/H2_Coulomb.'+str(iR), 'w')
	print >> fo, 'Uc=zeros(('+str(len(Sol))+','+str(len(Sol))+','+str(len(Sol))+','+str(len(Sol))+'),dtype=float)'
	for i,(i1,i2) in enumerate(ist):
	    (R1,Ene1,m1,p1,A1) = Sol[i1]
	    (R2,Ene2,m2,p2,A2) = Sol[i2]
	    for j,(i3,i4) in enumerate(ist):
	        (R3,Ene3,m3,p3,A3) = Sol[i3]
	        (R4,Ene4,m4,p4,A4) = Sol[i4]
	        if abs(Uc[i1,i2,i3,i4])>1e-7:
	            print >> fo, 'Uc[', i1, ',', i2, ',', i3, ',', i4, '] =', Uc[i1,i2,i3,i4]
	
