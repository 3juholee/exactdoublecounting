import sys
import os
from scipy import *
from scipy import optimize
from scipy import special
from scipy import integrate
from pylab import *
import BabHil
import poly as polyx
FORTRAN=True

def alp(l,m,p,A,R1):
    q = p
    l2 = l**2
    m2 = m**2
    l22 = 4*l2-1.
    q2 = q**2
    alp0 = (l2-m2)*(R1**2 - 4*q2*l2)/l22
    alp1 = A - q2 + l2-l
    alp3 = 8*p*l2 * (l2-m2)/l22
    return (alp0,alp1,alp3)

def gam(k,m,p,A,R2):
    p2 = p**2
    r2p = R2/(2*p)
    gam0 = (k+1)*(k+m+1)*(k+1-r2p)*(k+m+1-r2p)
    gam1 = A - p2 + R2 - (m+1)*(2*p+1-r2p) - 2*k*(k+m+2*p+1-r2p)
    gam2 = 2*p + (2*k+m+1)*(2+r2p/p)
    gam3 = (k+1)*(k+m+1)*(r2p/p)*(2*r2p-2*k-2-m)
    return (gam0,gam1,gam2,gam3)

def downBaberF(p,A,m,Kmax,R1,start1=1e-10):
    F2 = 0
    F1 = start1
    for k in range(Kmax,-1,-1):
        (alp0,alp1,alp3) = alp(m+k+1,m,p,A,R1)
        F0 = -alp1*F1 - alp0*F2
        F2 = F1
        F1 = F0
    return F0

def downBaberD(p,A,m,Kmax,R1,start1=1e-10):
    F2 = 0
    F1 = start1
    dFA2 = 0
    dFA1 = 0
    dFp2 = 0
    dFp1 = 0
    for k in range(Kmax,-1,-1):
        (alp0,alp1,alp3) = alp(m+k+1,m,p,A,R1)
        F0 = -alp1*F1 - alp0*F2
        dFA0 = -alp1*dFA1-alp0*dFA2-F1
        dFp0 = -alp1*dFp1-alp0*dFp2+2*p*F1+alp3*F2
        F2 = F1
        dFA2 = dFA1
        dFp2 = dFp1
        F1 = F0
        dFA1 = dFA0
        dFp1 = dFp0
    return (F0, dFA0, dFp0)


def downHyllerF(p,A,m,Kmax,R2,start1=1e-10):
    G2 = 0
    G1 = start1
    for k in range(Kmax,-1,-1):
        (gam0,gam1,gam2,gam3) = gam(k,m,p,A,R2)
        G0   = -gam1*G1   -gam0*G2
        G2 = G1
        G1 = G0
    return G0


def downHyllerD(p,A,m,Kmax,R2,start1=1e-10):
    G2 = 0
    G1 = start1
    dGA2 = 0
    dGA1 = 0
    dGp2 = 0
    dGp1 = 0
    for k in range(Kmax,-1,-1):
        (gam0,gam1,gam2,gam3) = gam(k,m,p,A,R2)
        G0   = -gam1*G1   -gam0*G2
        dGA0 = -gam1*dGA1 -gam0*dGA2 - G1
        dGp0 = -gam1*dGp1 -gam0*dGp2 + gam2*G1 + gam3*G2
        G2 = G1
        dGA2 = dGA1
        dGp2 = dGp1
        G1 = G0
        dGA1 = dGA0
        dGp1 = dGp0
    return (G0, dGA0, dGp0)


def Babers(p,A,m,Kmax,R1,start1=1e-10):
    def alphas(l,m,p,A,R1):
        q = p
        alp0 = (l-m)*(R1-2*q*l)/(2*l-1)
        alp1 = A - q**2 + l*(l+1)
        alp2 = (l+m+1)*(2*q*(l+1)+R1)/(2*l+3)
        return (alp0,alp1,alp2)
   
    a2 = 0
    a1 = start1
    aks=[]
    for k in range(Kmax,0,-1):
        (alp0,alp1,alp2) = alphas(m+k,m,p,A,R1)
        if alp0 == 0: 
         print 'q =', p
         print 'R1 =', R1
         print 'k =', k
         sys.exit()
        a0 = -a1*alp1/alp0 - a2*alp2/alp0
        a2 = a1
        a1 = a0
        aks.append(a0)
    aks.reverse()
    aks = array(aks)/aks[0]
    return aks

def Hyllers(p,A,m,Kmax,R2,start1=1e-10):
    def gammas(k,m,p,A,R2):
        p2 = p**2
        r2p = R2/(2*p)
        gam0 = k*(k+m-r2p)
        gam1 = A - p2 + R2 - (m+1)*(2*p+1-r2p) - 2*k*(k+m+2*p+1-r2p)
        gam2 = (k+m+1)*(k+1-r2p)
        return (gam0,gam1,gam2)
    
    C2 = 0
    C1 = start1
    Cks=[]
    for k in range(Kmax,0,-1):
        (gam0,gam1,gam2) = gammas(k,m,p,A,R2)
        C0 = -C1*gam1/gam0  - C2*gam2/gam0
        C2 = C1
        C1 = C0
        Cks.append(C0)
    Cks.reverse()
    Cks = array(Cks)/Cks[0]
    return Cks

def BaberTest(p,A,m,Kmax,R1,start1=1e-10):
    F2 = 0
    F1 = start1
    Fks=[]
    for k in range(Kmax,-1,-1):
        (alp0,alp1,alp3) = alp(m+k+1,m,p,A,R1)
        F0 = -alp1*F1 - alp0*F2
        Fks.append(F1)
        F2 = F1
        F1 = F0
    Fks.reverse()
    q=p
    aks=[1.0]
    prod=1
    for k in range(1,Kmax):
        prod *= k*(R1-2*q*(m+k))/(2*(m+k)-1)
        aks.append( Fks[k]*prod/Fks[0] )
    return array(aks)

def HyllerTest(p,A,m,Kmax,R2,start1=1e-10):
    G2 = 0
    G1 = start1
    Gks=[]
    for k in range(Kmax,-1,-1):
        (gam0,gam1,gam2,gam3) = gam(k,m,p,A,R2)
        G0   = -gam1*G1   -gam0*G2
        Gks.append(G1)
        G2 = G1
        G1 = G0
    Gks.reverse()
    
    Cks=[1.0]
    prod=1
    for k in range(1,Kmax):
        prod *= k*(k+m-R2/(2*p))
        Cks.append( Gks[k]*prod/Gks[0] )
    return array(Cks)

class BaberHyller:
    def __init__(self,m,Kmax,R,Za,Zb):
        self.m = m
        self.R = R
        self.R1 = R*(Za-Zb)
        self.R2 = R*(Za+Zb)
        self.Kmax = Kmax
        self.start1 = 1e-10
    def F(self, x):
        p = x[0]
        A = x[1]
        if (FORTRAN):
            F = BabHil.baberf(p,A,self.m,self.Kmax,self.R1)
            G = BabHil.hyllerf(p,A,self.m,self.Kmax,self.R2)
        else:
            F = downBaberF(p,A,self.m,self.Kmax,self.R1)
            G = downHyllerF(p,A,self.m,self.Kmax,self.R2)
        return (F, G)
    def D(self, x):
        p = x[0]
        A = x[1]
        if (FORTRAN):
            (F, dFA, dFp) = BabHil.baberd(p,A,self.m,self.Kmax,self.R1)
            (G, dGA, dGp) = BabHil.hyllerd(p,A,self.m,self.Kmax,self.R2)
        else:
            (F, dFA, dFp) = downBaberD(p,A,self.m,self.Kmax,self.R1)
            (G, dGA, dGp) = downHyllerD(p,A,self.m,self.Kmax,self.R2)
        return array([[dFp,dFA],[dGp,dGA]])

    def StartGuess(self, n,l,m):
        def hx(l,m,p):
            if l==0: return 0
            else: return (l**2-m**2)*((self.R1/(2*p))**2-l**2)/(4*l**2-1)/l
            
        p=self.R2/(2.*n)
        A=-l*(l+1)+2*p**2*(hx(l+1,m,p)-hx(l,m,p)+0.5)
        return (p,A)


def FindSolution(n,l,m, Kmax,Za,Zb):
    R_start = 1.2
    R_min = 0.6
    R_max = 12
    N_tot = 500
    N1 = int((R_max-R_start)/(R_max-R_min)*N_tot)
    N2 = int((R_start-R_min)/(R_max-R_min)*N_tot)
    sol1=[]
    for iR in range(N1):
        R = R_start + iR/(N1-1.)*(R_max-R_start)
        BH = BaberHyller(m,Kmax,R,Za,Zb)
        if iR==0:  (p,A) = (R*(Za+Zb)/(2.*n), -l*(l+1))   # He-atom for small R!
        (p,A)= optimize.fsolve(BH.F, [p,A], fprime=BH.D, full_output=0)
        Ene = -p**2*4/R**2 + 2./R
        sol1.append([R, Ene, p, A])
    sol2=[]
    for iR in range(N2):
        R = R_start - (iR+1.)/N2*(R_start-R_min)
        BH = BaberHyller(m,Kmax,R,Za,Zb)
        if iR==0:  (p,A) = (R*(Za+Zb)/(2.*n), -l*(l+1))   # He-atom for small R!
        (p,A)= optimize.fsolve(BH.F, [p,A], fprime=BH.D, full_output=0)
        Ene = -p**2*4/R**2 + 2./R
        if iR>0 : sol2.append([R, Ene, p, A])
    sols = sol2[::-1]+sol1[:]
    return array(sols)
    

def NormBaber(aks,bks,m,q):
    "Computes the normalization constants N and N' of Baber function M!"
    dsum1=0.0
    for l in range(m,len(aks)+m):
        norm_Plm = 2./(2.*l+1)*special.gamma(l+m+1)/special.gamma(l-m+1)
        dsum1 += aks[l-m]*bks[l-m]*(-1)**l*norm_Plm

    x0=0.5
    dsuma_0 = sum([aks[l-m]*special.lpmv(m,l,x0) for l in range(m,len(aks)+m)])
    dsumb_0 = sum([bks[l-m]*(-1)**l*special.lpmv(m,l,x0) for l in range(m,len(aks)+m)])
    Np_over_N = (-1)**m*exp(-2*x0*q)*dsuma_0/dsumb_0
    Norma = sqrt(abs(1./(Np_over_N * dsum1)))
    Normb = Np_over_N * Norma
    #print 'Np_over_N=', Np_over_N, 'dsuma=', dsuma_0, 'dsumb=', dsumb_0, 'dsum1=', dsum1
    
    return (Norma,Normb)

def M2mu2(aks,bks,m):
    "Computes the integral  Integrate[mu^2 * M(mu)^2, {mu,-1,1}] "
    def l1_e_l2(l,m):
        rs = ((l-m+1.)/(2*l+1.))**2*2/(2*l+3.)*special.gamma(l+m+2.)/special.gamma(l-m+2.)
        if l>m: rs+=((l+m)/(2*l+1.))**2*2/(2*l-1.)*special.gamma(l+m)/special.gamma(l-m)
        return rs
    def l1_2_l2(l,m):
        return 2./((2*l-1.)*(2*l+1.)*(2*l+3.))*special.gamma(l+m+2)/special.gamma(l-m)
    
    dsum1=0.0
    for l in range(m,len(aks)+m):
        dsum1 += (-1)**l*aks[l-m]*bks[l-m]*l1_e_l2(l,m)
    dsum2=0.0
    for l in range(m+1,len(aks)+m-1):
        dsum2 += (-1)**l*(aks[l-m-1]*bks[l-m+1]+aks[l-m+1]*bks[l-m-1])*l1_2_l2(l,m)
    return (-1)**m*(dsum1-dsum2)


def testM(x,m,q,aks):
    Mx = 0.0
    for l in range(m,len(aks)+m):
        Mx += aks[l-m]*special.lpmv(m,l,x)
    return x**2*(exp(-q*x) * Mx)**2


def NormHyller(x,m,p,Cks,M2mu2):
    #csum = sum([Cks[k]*special.eval_genlaguerre(m+k,m,x) for k in range(len(Cks))]) # original
    csum = sum([Cks[k]*special.eval_genlaguerre(k,m,x) for k in range(len(Cks))])
    lambd = 1+x/(2*p)
    return (lambd**2-1)**m*exp(-x)*csum**2*(lambd**2-M2mu2)/(2*p)

def NormHyller_dirty(x,m,p,Cks,M2mu2):
    alpha = m+0.0
    #print 'alpha',alpha
    max_n = m+len(Cks)-1
    Lag = polyx.Laguerre_dirty(x,alpha,max_n)
    #csum = sum([Cks[k]*Lag[m+k] for k in range(len(Cks))]) # original
    csum = sum([Cks[k]*Lag[k] for k in range(len(Cks))])
    lambd = 1+x/(2*p)
    return (lambd**2-1)**m*exp(-x)*csum**2*(lambd**2-M2mu2)/(2*p)

def GiveExpansion((R,E,p,A),m,Kmax):
	R1 = R*(Za-Zb)
	R2 = R*(Za+Zb)
	#print 'p=', p, 'A=', A
	
	#### Finding coefficients for M(mu) due to Baber
	aks = Babers(p,A,m,Kmax,R1)
	bks = Babers(p,A,m,Kmax,-R1)
	(Norma,Normb) = NormBaber(aks,bks,m,p)
	aks = aks[:]*Norma
	bks = bks[:]*Normb
	m2mu2 = M2mu2(aks,bks,m)
	#### Finding coefficients for Lambda(lambda)
	Cks = Hyllers(p,A,m,Kmax,R2)
	
	
	Romb=False
	if (Romb):
		Nx = 2**9+1 
		Xx = linspace(0,30,Nx)
		Tt = exp(Xx)
		Hyller_array = zeros((len(Tt)),float)
		for it,t in enumerate(Tt):
			Hyller_array[it] = NormHyller(t,m,p,Cks,m2mu2)
		Hy_Tt = Hyller_array*Tt
		NormH = integrate.romb(Hy_Tt,Xx,Xx[1]-Xx[0])

	else:
	    NormH = integrate.quad(NormHyller_dirty,0,Inf,args=(m,p,Cks,m2mu2))[0]
	    
	#print 'NormH=', NormH
	#print 'NormH=', NormH
	
	Cks = Cks[:]/sqrt(NormH)
	return (aks, Cks, NormH)
    

if __name__ == '__main__':
    Kmax=20
    Za = 1
    Zb = 1
    
#enter the input(initial states)
    execfile('input.py')
			
    Solutions={}
    for (n,l,m) in n_l_m:
        sols = FindSolution(n,l,abs(m),Kmax,Za,Zb)
        Solutions[(n,l,m)]=sols
        print 'Done n=',n,'l=',l,'m=',m
    
    Rs = Solutions[n_l_m[0]][:,0]
		# 0th column of sol : range R

    if not os.path.exists('H2+'): os.makedirs('H2+')
    
    for nkey in n_l_m:
        (n,l,m) = nkey
        sols = Solutions[nkey]
        fi = open('H2+/H2_n_'+str(n)+'_l_'+str(l)+'_m_'+str(m), 'w')
        print >> fi, '# n=', n, 'l=', l, 'm=', m
        print >> fi, '#  R  Energy   p    A'
        for i,s in enumerate(sols):
            print >> fi, "%20.12g %20.12g %20.12g %20.12g" % tuple(sols[i])
    
    
    execfile('which.py')
    which.sort()

    for iR in which:
        fr = open('H2+/H2_solution.'+str(iR),'w')
        print >> fr, '# sol=[R,Ene,m,p,A]  aks=[...] Cks=[...]'
        norms=[]
        for nkey in n_l_m:
            (n,l,m)=nkey
            sol = Solutions[nkey][iR]
            (aks, Cks, NormH) = GiveExpansion(sol,abs(m),Kmax)
            norms.append(NormH)
            [R, Ene, p, A] = sol
            print >> fr, 'sol=', [R,Ene,m,p,A], '; aks=', aks.tolist(), '; Cks=', Cks.tolist()
        print 'iR=', iR, 'R=', Rs[iR]#, norms
