#popen3 -> subprocess.Popen for python ver 2.6

import subprocess
import sys
import os
import shutil
import re
from scipy import *
from scipy import weave
from scipy import linalg
from scipy import optimize
from scipy import interpolate
from h2utils import *
from cix import *
import copy
import time

def create_log_mesh(sigdata, nom, ntail_):
    """Creates logarithmic mesh on Matsubara axis
       Takes first istart points from mesh om and
       the rest of om mesh is replaced by ntail poinst
       redistribued logarithmically.
       Input:
           om      -- original long mesh
           istart  -- first istart points unchanged
           ntail   -- tail replaced by ntail points only
       Output:
           som             -- smaller mesh created from big om mesh
           sSig[nom,nc]    -- Sig on small mesh
       Also computed but not returned:
           ind_om  -- index array which conatins index to
                      kept Matsubara points
    """
    om = sigdata[0]
    
    istart = min(nom, len(om))
    ntail = min(ntail_, len(om)-istart)
        
    istart = min(nom,len(om))
    ntail = min(ntail, len(om)-istart)

    ind_om=[]
    alpha = log((len(om)-1.)/istart)/(ntail-1.)
    for i in range(istart):
        ind_om.append(i)
    for i in range(ntail):
        t = int(istart*exp(alpha*i)+0.5)
        if (t != ind_om[-1]):
            ind_om.append(t)

    ind_oms_equal = [[0]]
    for it in range(1,len(ind_om)-1):
        istart = int(0.5*(ind_om[it-1]+ind_om[it])+0.51)
        iend = int(0.5*(ind_om[it]+ind_om[it+1])-0.01)
        equal = [i for i in range(istart,iend+1)]
        ind_oms_equal.append(equal)
    istart = int(0.5*(ind_om[-2]+ind_om[-1])+0.51)
    equal = [i for i in range(istart,ind_om[-1]+1)]
    ind_oms_equal.append(equal)

    oms_equal=[]
    for ind in ind_oms_equal:
        oms_equal.append( array([om[i] for i in ind]) )
    #print ind_oms_equal
    #print oms_equal
    
    ssigdata = zeros( (shape(sigdata)[0], len(ind_om)), dtype=float )
    for i in range(len(ind_om)):
        ssigdata[:,i] = sigdata[:,ind_om[i]]

    return (ssigdata,oms_equal,ind_om)


def ferm(x):
    if x>300: return 0.0
    if x<-300: return 1.0
    return 1/(exp(x)+1.0)



def MatsubaraSum(mu,Ew,om,beta):
    sm=0.0
    for iw,w in enumerate(om): sm += real(1/(w*1j+mu-Ew[iw]))
    sm = 2.*sm/beta+0.5
    # Sommerfeld expansion to correct for finite sum
    # We need to add sum_{iw>Omega0} 1/(iw-E0)
    # which is Integrate[f(x) * Omega0/(Omega0**2+(x-E0)**2)/pi,{x,-oo,oo}]
    E0 = Ew[-1].real
    Omega0 = om[-1]+pi/beta
    sm -= 1./pi*arctan2(E0,Omega0)
    sm += pi/(3*beta**2)*Omega0*E0/(Omega0**2+E0**2)**2
    sm -= 7*pi**3/(15*beta**4)*Omega0*E0*(Omega0**2-E0**2)/(Omega0**2+E0**2)**4
    return sm


def LogMatsubaraSum(mu,Ew,oms,oms_equal,beta):
    code="""
       #include <complex>
       using namespace std;
       double sm=0.0;
       double cmu=mu;
       complex<double> i(0.0,1.0);
       for (int iband=0; iband<Ec.size(); iband++){
           complex<double> cEc = Ec(iband);
           for (int j=0; j<equal.size(); j++){
               double w = equal(j);
               complex<double> G = 1.0/(i*w+cmu-cEc);
               sm += G.real();
           }
       }
       return_val = sm;
    """
    PYTHON=False
    if (PYTHON):
        sm=0.0
        for iw in range(len(oms)):
            for Eb in Ew[iw]:
                for w in oms_equal[iw]:
                    sm += real(1/(w*1j+mu-Eb))
    else:
        sm=0.0
        for iw in range(len(oms)):
            equal = oms_equal[iw]
            Ec = Ew[iw]
            sm += weave.inline(code, ['equal', 'mu', 'Ec'], type_converters=weave.converters.blitz, compiler = 'gcc')
        
    sm *= 2./beta
    for Eb in Ew[-1]:
        # Sommerfeld expansion to correct for finite sum
        # We need to add sum_{iw>Omega0} 1/(iw-E0)
        # which is Integrate[f(x) * Omega0/(Omega0**2+(x-E0)**2)/pi,{x,-oo,oo}]
        E0 = Eb.real-mu
        Omega0 = oms[-1]+pi/beta
        sm += 0.5  # because of 1/omega
        sm -= 1./pi*arctan2(E0,Omega0)
        sm += pi/(3*beta**2)*Omega0*E0/(Omega0**2+E0**2)**2
        sm -= 7*pi**3/(15*beta**4)*Omega0*E0*(Omega0**2-E0**2)/(Omega0**2+E0**2)**4
    return sm


def LogMatsubaraSum2(mu,Ew,oms,oms_equal,beta):
    code="""
       #include <complex>
       using namespace std;
       double sm=0.0;
       double cmu=mu;
       complex<double> i(0.0,1.0);
       for (int iband=0; iband<Ec.size(); iband++){
           complex<double> cEc = Ec(iband);
           for (int j=0; j<equal.size(); j++){
               double w = equal(j);
               complex<double> G = 1.0/(i*w+cmu-cEc);
               sm += G.real();
           }
       }
       return_val = sm;
    """
    PYTHON=True
    if (PYTHON):
        sm=zeros(shape(Ew)[1],dtype=float)
        for iw in range(len(oms)):
            for iE,Eb in enumerate(Ew[iw]):
                #print 'Eb=', Eb, 'mu=', mu, 'beta=', beta, oms_equal[iw]
                for w in oms_equal[iw]:
                    sm[iE] += real(1/(w*1j+mu-Eb))
                    #if iE==0: print w, Eb, mu, real(1/(w*1j+mu-Eb))
        #print 'sm0=', sm
    else:
        sm=0.0
        for iw in range(len(oms)):
            equal = oms_equal[iw]
            Ec = Ew[iw]
            sm += weave.inline(code, ['equal', 'mu', 'Ec'], type_converters=weave.converters.blitz, compiler = 'gcc')
        
    sm *= 2./beta
    for iE,Eb in enumerate(Ew[-1]):
        # Sommerfeld expansion to correct for finite sum
        # We need to add sum_{iw>Omega0} 1/(iw-E0)
        # which is Integrate[f(x) * Omega0/(Omega0**2+(x-E0)**2)/pi,{x,-oo,oo}]
        E0 = Eb.real-mu
        Omega0 = oms[-1]+pi/beta
        sm[iE] += 0.5  # because of 1/omega
        #if iE==0: print 'sm1=', sm[iE]
        sm[iE] -= 1./pi*arctan2(E0,Omega0)
        sm[iE] += pi/(3*beta**2)*Omega0*E0/(Omega0**2+E0**2)**2
        sm[iE] -= 7*pi**3/(15*beta**4)*Omega0*E0*(Omega0**2-E0**2)/(Omega0**2+E0**2)**4
        #if iE==0: print 'sm2=', sm[iE]
    return sm



def cmp_UC(Uc,T):
    # Transforming U of 2x2 system to local basis
    UC = zeros((2,2,2,2),dtype=float)
    for i1 in range(2):
        for i2 in range(2):
            for i3 in range(2):
                for i4 in range(2):
                    for a in range(2):
                        for b in range(2):
                            for c in range(2):
                                for d in range(2):
                                    UC[i1,i2,i3,i4] += T[a,i1]*T[b,i2]*T[c,i3]*T[d,i4]*Uc[a,b,c,d]
    return UC



def Cmp_Gloc_Delta(lsigdata,s_oo,mu,Edc,Ener,TProj):
    N = len(Ener)
    Eimp = sum([TProj[0,i]**2 * (Ener[i]-mu) for i in range(N)])-Edc

    Glocal=zeros((len(lsigdata[0]),3),dtype=float)
    Delta=zeros((len(lsigdata[0]),3),dtype=float)
    Glocal[:,0] = lsigdata[0]
    Delta[:,0] = lsigdata[0]
    for iw,w in enumerate(lsigdata[0]):
        iomega = w*1j
        #iomega = w
        Sigm = (lsigdata[1,iw]+lsigdata[2,iw]*1j)+s_oo[0]
        #print 'shape(lsigdata)', shape(lsigdata), 'shape(lsigdata[1,iw])=', lsigdata[1,iw], lsigdata[2,iw], shape(s_oo)
        Sig = identity(2)*(Sigm-Edc)
        Ginv = zeros((N,N),dtype=complex)                
        for i in range(N): Ginv[i,i] = iomega+mu-Ener[i]
        Ginv -= TProj.T * Sig * TProj
        Gc = linalg.inv(Ginv)
        Gloc = (TProj * Gc * TProj.T)[0,0]
        #print 'shape(Gloc)=',shape(Gloc)
        #print 'shape(Gc)=', shape(Gc)
        Delt = iomega-Eimp-Sigm-1/Gloc
        #print 'shape(Delt)=', shape(Delt)
        #print 'shape(iomega)=', iomega
        #print 'shape(Eimp)=', Eimp
        #print 'shape(Sigm)=', Sigm
        
        Glocal[iw,1:] = array([Gloc.real,Gloc.imag])

        #print 'Delta=', Delt.real, Delt.imag
        #print array([Delt.real,Delt.imag])
        #print 'shape(Delta)=', shape(Delta[iw,1:])
        
        Delta [iw,1:] = array([Delt.real,Delt.imag])
    return (Eimp,Glocal,Delta)


def Cmp_eks(lsigdata,s_oo,Edc,Ener,TProj):

    N = len(Ener)
    eks=zeros((len(lsigdata[0]),N),dtype=complex)
    for iw,w in enumerate(lsigdata[0]):
        iomega = w*1j
        Sigm = (lsigdata[1,iw]+lsigdata[2,iw]*1j)+s_oo[0]
        Sig = identity(2)*(Sigm-Edc)
        ES = zeros((N,N),dtype=complex)
        for i in range(N): ES[i,i] = Ener[i]
        ES += TProj.T * Sig * TProj
        eks[iw] = linalg.eigvals(ES)
    return eks

def Cmp_eVks(lsigdata,s_oo,Edc,Ener,TProj):
    
    N = len(Ener)
    eks=zeros((len(lsigdata[0]),N),dtype=complex)
    Ars=zeros((len(lsigdata[0]),N,N),dtype=complex)
    for iw,w in enumerate(lsigdata[0]):
        iomega = w*1j
        if iw!=N-1:
            Sigm = (lsigdata[1,iw]+lsigdata[2,iw]*1j)+s_oo[0]
        else:
            Sigm = lsigdata[1,iw]+s_oo[0]
            
        Sig = identity(2)*(Sigm-Edc)
        ES = zeros((N,N),dtype=complex)
        for i in range(N): ES[i,i] = Ener[i]
        ES += TProj.T * Sig * TProj

        if iw!=N-1:
            (es, Ar) = linalg.eig(ES)
        else:
            (es, Ar) = linalg.eigh(ES)
        
        eks[iw] = es
        Ars[iw,:,:] = Ar
    return (eks,Ars)
    
def Cmp_mu(N0,a,b,eks,oms,oms_equal,beta):

    def Density(mu,N0,eks,oms,oms_equal,beta):
        return LogMatsubaraSum(mu,eks,oms,oms_equal,beta)-N0

    return optimize.brentq(Density,a,b,args=(N0,eks,oms,oms_equal,beta))
    

def Cmp_DensityMatrix(eks,Ars,mu,oms,oms_equal,beta):

    code1="""
      using namespace std;
      complex<double> i(0,1);
      complex<double> sm=0.0;
      for (int jw=0; jw<equal.size(); jw++){
         double w = equal(jw);
         sm += 1.0/(w*i+mu-Eb);
      }
      return_val = sm;
    """
    code2="""
      using namespace std;
      for (int i1=0; i1<N; i1++){
         for (int i2=0; i2<N; i2++){
            for (int l=0; l<N; l++){
               Glc(i1,i2) += Ar(i1,l)*gii(l)*Al(l,i2);
            }
         }
      }
    """
    N = len(eks[0])
    Nc0 = zeros((N,N),dtype=complex)
    for iw in range(len(oms)):
        Ar = Ars[iw]
        Al = linalg.inv(Ar)
        gii = zeros((N,),dtype=complex)
        
        for l in range(N):
            Eb = eks[iw,l]
            equal = oms_equal[iw]
            sm = weave.inline(code1, ['equal','mu','Eb'], type_converters=weave.converters.blitz, compiler = 'gcc')
            
            if iw==len(oms)-1:
                # Sommerfeld expansion to correct for finite sum
                # We need to add sum_{iw>Omega0} 1/(iw-E0)
                # which is Integrate[f(x) * Omega0/(Omega0**2+(x-E0)**2)/pi,{x,-oo,oo}]
                E0 = Eb.real-mu
                Omega0 = oms[-1]+pi/beta
                sm_oo = 0.5  # because of 1/omega
                sm_oo -= 1./pi*arctan2(E0,Omega0)
                sm_oo += pi/(3*beta**2)*Omega0*E0/(Omega0**2+E0**2)**2
                sm_oo -= 7*pi**3/(15*beta**4)*Omega0*E0*(Omega0**2-E0**2)/(Omega0**2+E0**2)**4
                sm += 0.5*beta*sm_oo
            gii[l] = sm
    
        Glc = zeros((N,N),dtype=complex)
        weave.inline(code2, ['Glc','N','Ar','Al','gii'], type_converters=weave.converters.blitz, compiler = 'gcc')
        Nc0 += (Glc+transpose(conjugate(Glc)))/beta
        
    return Nc0
        



class IMP:
    def __init__(self,nom,beta,RUUN=True,wbroad=0.03,kbroad=0.15):
        self.checkCausality = False
        self.broad = True
        
        self.params = {"exe":   ["./ctqmc",          "# Path to executable"],
                       "workdir":["imp",             "# Working directory for impurity"],
                       "Delta": ["Delta.dat",        "# Input bath function hybridization"],
                       "Sig":   ["Sig.out",          "# Output self-energy"],
                       "cix":   ["one_band.imp",     "# Input file with atomic state"],
                       "U":     [0,                  "# Coulomb repulsion (F0)"],
                       "mu":    [0,                  "# Chemical potential"],
                       "beta":  [100,                "# Inverse temperature"],
                       "M" :    [30e6,               "# Number of Monte Carlo steps"],
                       "nom":   [150,                "# number of Matsubara frequency points to sample"],
                       "nomD":  [50,                 "# number of Matsubara points using the Dyson Equation"],
                       "aom":   [5,                  "# number of frequency points to determin high frequency tail"],
                       "tsample":[300,               "# how often to record the measurements" ],
                       "maxNoise":[1e100,            "# maximum allowed noise is large in simple run"]}
        self.wbroad=wbroad                     # Broadening of the hybridization function
        self.kbroad=kbroad                     # Broadening of the hybridization function

        self.dir = self.params['workdir'][0]
        self.q_exe = self.params['exe'][0]
        self.Signame = self.params['Sig'][0]
        self.PARAMS = 'PARAMS'
        self.fh_info = sys.stdout
        self.mpi_prefix = ''
        mpifile = 'mpi_prefix.dat'
        if os.path.isfile(mpifile):
            self.mpi_prefix = open(mpifile, 'r').next().strip()
            print "DmftEnvironment: mpi_prefix.dat exists -- running in parallel mode."
            print "  ", self.mpi_prefix
        else:
            print "DmftEnvironment: mpi_prefix.dat does not exist -- running in single-processor mode."
        
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        
        # creating impurity cix file
        f = open(self.dir+'/'+self.params['cix'][0], 'w')
        print >> f, one_band_cix
        f.close()

        root = os.getenv('WIEN_DMFT_ROOT')
        shutil.copy2(root+'/ctqmc', self.dir)
        shutil.copy2(root+'/broad', self.dir)
        
        self.params['nom'][0]=nom
        self.params['beta'][0]=beta
        
    def Run(self,U,mu_QMC,omega_large,ind_om,Glocal,Delta,s_oo):

        savetxt(self.dir+'/Gloc.dat',Glocal)
        
        # Interpolating Delta on entire Matsubara mesh
        Dreal = interpolate.UnivariateSpline(Delta[:,0],Delta[:,1],s=0)
        Dimag = interpolate.UnivariateSpline(Delta[:,0],Delta[:,2],s=0)
        Dr = [Dreal(x) for x in omega_large]
        Di = [Dimag(x) for x in omega_large]
        Delta2 = transpose(array([omega_large,Dr,Di]))
        savetxt(self.dir+'/Delta.dat',Delta2)
        
        # Creates input file (PARAMS) for CT-QMC solver
        self.params['U'][0]=U
        self.params['mu'][0]=mu_QMC
        f = open(self.dir+'/PARAMS', 'w')
        print >> f, '# Input file for continuous time quantum Monte Carlo'
        for p in self.params.keys():
            print >> f, p, self.params[p][0], '\t', self.params[p][1]
        f.close()
       
        # Below we execute ctqmc
        cmd = 'cd '+self.dir+'; '+self.mpi_prefix+' '+self.q_exe+' '+self.PARAMS+' > nohup_imp.out 2>&1 '
        print cmd
        print 'Running ---- ctqmc -----'

        RUUN = True
        if (RUUN):
            p = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE, close_fds=True)
            stdin, stdout, stderr = (p.stdin,p.stdout,p.stderr)
            print >> self.fh_info, stdout.read(), stderr.read()

        
        # Copying to Sig.outb
        Signame = self.dir+'/'+self.Signame
        shutil.copy2(Signame, Signame+'b')

        nf=0.; TrSigmaG=0.
        first_line = open(Signame,'r').readline()
        m=re.search('nf=(\d*\.\d*)',first_line)
        if m is not None:
            nf = float(m.group(1))
        m=re.search('TrSigmaG=(\d*\.\d*)',first_line)
        if m is not None:
            TrSigmaG = float(m.group(1))


        Sig = loadtxt(Signame)
        if (RUUN):
            s_oo = copy.deepcopy(Sig[-1,1::2])
            for l in range(len(s_oo)): Sig[:,2*l+1]-=s_oo[l]
            
        if self.broad:
            
            if self.checkCausality:
                for i in range(shape(Sig)[1]/2):
                    for iw in range(len(Sig)):
                        if Sig[iw,2+2*i]>0: Sig[iw,2+2*i]=0.0
            
            savetxt(Signame+'w',Sig)
            
            # adds the first line for header at the beginning
            f0=open(Signame+'w','r')
            dat=f0.readlines()
            f0.close()
            f1=open(Signame+'t','w')
            f1.writelines([first_line]+dat)
            f1.close()
            shutil.move(Signame+'t',Signame+'w')
            
            # Broadening the output self-energy to reduce qmc noise
            cmd = 'cd '+self.dir+'; ./broad -w '+str(self.wbroad)+' -k '+str(self.kbroad)+' '+self.Signame+'w'+' > '+self.Signame
            print cmd
            
            if (RUUN):
                p = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE, close_fds=True)
                stdin, stdout, stderr = (p.stdin,p.stdout,p.stderr)
                print >> self.fh_info, stdout.read(), stderr.read()
            
            Sig = loadtxt(Signame)

            
        Sigs = zeros( (len(ind_om), shape(Sig)[1]), dtype=float )
        for i in range(len(ind_om)):
            Sigs[i,:] = Sig[ind_om[i],:]
        
        return (transpose(Sigs), s_oo, nf, TrSigmaG)
    
    def copyfiles(self,itt):
        shutil.copy2(self.dir+'/'+self.params['Sig'][0], self.dir+'/'+self.params['Sig'][0]+'.'+str(itt))
        shutil.copy2(self.dir+'/nohup_imp.out', self.dir+'/nohup_imp.out.'+str(itt))
        shutil.copy2(self.dir+'/Gf.out', self.dir+'/Gf.out.'+str(itt))
        shutil.copy2(self.dir+'/Gloc.dat', self.dir+'/Gloc.dat.'+str(itt))
        

    
def cmp_Gloc(iw,oms,Sigs,s_oo,Edc,mu,Ener,TProj):
    iomega = oms[iw]*1j
    Sigm = (Sigs[1,iw]+Sigs[2,iw]*1j)+s_oo[0]
    Sig = identity(2)*(Sigm-Edc)
    ES = zeros((N,N),dtype=complex)
    for i in range(N): ES[i,i] = Ener[i]
    ES += TProj.T * Sig * TProj
    Gloc = linalg.inv( (iomega+mu)*identity(N)-ES )
    return Gloc
    
if __name__ == '__main__':
    
    if len(sys.argv)<=1:
        print 'Give input integer iR!'
        sys.exit(1)
    
    DMFT=True
    RUUN=True

    Etol=1e-6
    Nomega = 5000
    beta = 100
    nom = 150
    ntail = 100
    gamma = 0.0
    Noccup=2.0  # occupancy for H2
    mu_a, mu_b = -2., 2.
    mixr = 1.0
    Nitt = 100 # maximum number of iterations
    Nitt = 1 # maximum number of iterations
    
    omega_large = linspace(pi/beta,pi*(2*Nomega-1)/beta,Nomega)
    print omega_large[0], omega_large[-1]
    sigdata = zeros((3,len(omega_large)),dtype=float)
    sigdata[0,:] = omega_large

    ctqmc = IMP(nom,beta)
    
    iR=int(sys.argv[1])
    print 'iR=', iR, 'DMFT=', DMFT
    
    (Sol,Ck,ak) = ReadSolution(iR)
    execfile('H2+/H2_Coulomb.'+str(iR))
    
    index=[]
    for i,(R,Ene,m,p,A) in enumerate(Sol):
        for deg in range(2*m+1): index.append(i)
    
    print 'All states size=', len(index), 'index=', index
    
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
    
    T=matrix([[1/sqrt(2.),1/sqrt(2.)],[1/sqrt(2.),-1/sqrt(2.)]])
    if (DMFT):
        UC = cmp_UC(Uc,T)
        U_local=UC[0,0,0,0]
        print 'U_onsite=', U_local


    if os.path.isfile(ctqmc.dir+'/'+ctqmc.Signame):
        # use the old file, if exists
        sold = loadtxt(ctqmc.dir+'/'+ctqmc.Signame).transpose()
        for i in range(1,len(sold)):
            Sf = interpolate.UnivariateSpline(sold[0,:],sold[i,:],s=0)
            sigdata[i,:] = Sf(omega_large)
        (s2,Edc) = loadtxt(ctqmc.dir+'/'+'s_oo')
        s_oo=[s2]
    else:
        sigdata[2,:] = -gamma
        s_oo = [0.5*U_local]
        Edc = 0.5*U_local
    
    Sigs,oms_equal,ind_om = create_log_mesh(sigdata, nom, ntail)
    
    
    H0 = zeros((len(index),len(index)),dtype=float)
    V2 = zeros((len(index),len(index)),dtype=float)
    Ham = zeros((len(index),len(index)),dtype=float)
    for i,ii in enumerate(index):
        H0[i,i]=Sol[ii][1]-1./R # because both electrons from H2+ problem include +2/R

    # Hartree and Fock potential
    UH = zeros((len(ist),len(ist)), dtype=float)
    UF = zeros((len(ist),len(ist)), dtype=float)
    for i, (i1,i2) in enumerate(ist):
        for j, (i3,i4) in enumerate(ist):
            UH[i,j] = Uc[i1,i3,i4,i2]
            UF[i,j] = Uc[i1,i3,i2,i4]
    
    Nc = zeros(len(ist),dtype=float)
    Nc0 = zeros(len(ist),dtype=float)
    VHF = zeros(len(ist),dtype=float)
    VH = zeros(len(ist),dtype=float)
    VF = zeros(len(ist),dtype=float)
    
    if os.path.exists('Nc.dat'):
        Nc = loadtxt('Nc.dat')
    else:
        # more appropriate for localized regime
        Nc[in12_1[(0,0)]]=1.0
        Nc[in12_1[(1,1)]]=1.0
   
    p_Ener=0
    
    for itt in range(Nitt):
        # Hartree-Fock
        VHF = dot(UH, Nc)-0.5*dot(UF, Nc)
        VH = dot(UH, Nc)
        VF = -0.5*dot(UF, Nc)
        for i, (i1,i2) in enumerate(in12):  V2[i1,i2]=VHF[i]
        Ham = H0+V2
        (Ener,Evec) = linalg.eigh(Ham)
        
        if (DMFT):
            
            Proj = Evec[0:2,:] # Projector
            print 'from-H^+-to-HF-Proj[0]=', "%8.5f "*(len(Evec)) % tuple(Proj[0,:])
            print 'from-H^+-to-HF-Proj[1]=', "%8.5f "*(len(Evec)) % tuple(Proj[1,:])
            Proj = matrix(Proj)
            TProj = T*Proj

            oms = Sigs[0]
            # computes Hartree-Fock energy levels
            eks=Cmp_eks(Sigs,s_oo,Edc,Ener,TProj)
            # computes the current chemical potential
            mu = Cmp_mu(Noccup/2.,mu_a,mu_b,eks,oms,oms_equal,beta)
            # computes Gloc & Delta
            (Eimp,Glocal,Delta) = Cmp_Gloc_Delta(Sigs,s_oo,mu,Edc,Ener,TProj)
            
            print 'Eimp=', Eimp, 'U=', U_local, 'Edc=', Edc, 's_oo=', s_oo
            print 'mu[0]=', mu
            #if itt==0: Eimp=-U_local/2.
            (Sigs, s_oo, n_imp, TrSigmaG) = ctqmc.Run(U_local,-Eimp,omega_large,ind_om,Glocal,Delta,s_oo)

            ctqmc.copyfiles(itt)
            
            Edc = U_local*n_imp/2.
            print 's_oo-Edc=', s_oo-Edc


            
            (eks,Ars) = Cmp_eVks(Sigs,s_oo,Edc,Ener,TProj)
            # recomputes the current chemical potential for charge neutrality
            mu = Cmp_mu(Noccup/2.,mu_a,mu_b,eks,oms,oms_equal,beta)
            print 'mu[1]=', mu

            # In hartree-fock basis
            Nt = 2.*Cmp_DensityMatrix(eks,Ars,mu,oms,oms_equal,beta)  # 2 due to spin
            
            Nt2 = matrix(Evec) * Nt * matrix(Evec).T  # In H2+ basis
            for i, (i1,i2) in enumerate(in12):  Nc0[i] = real(Nt2[i1,i2])
            
            TrSG = sum(Nc0*VHF) # evaluates Tr(Sigma*G) in H2+ basis
            EH = 0.5*sum(Nc0*VH)
            EF = 0.5*sum(Nc0*VF)
            Etot0 = sum([H0[i,i]*real(Nt2[i,i]) for i in range(len(H0))]) # In H2+ basis
            Etot = Etot0 + 0.5*TrSG # This is Etot2 = Tr(H0*G)+0.5*Tr(VHF*G)

            
            #print 'imag(Nc)=', sum(sum(abs(imag(Nt2))))
            
            Nc[:] = mixr*Nc0[:] + (1-mixr)*Nc[:]
            
            N_local = 0.5*(Nc[in12_1[(0,0)]]+Nc[in12_1[(0,1)]]+Nc[in12_1[(1,0)]]+Nc[in12_1[(1,1)]])

            Nf=n_imp
            #print 'Nf = ',Nf
            Edc = U_local*Nf/2.
            EHF_local = U_local*Nf**2/4. * 2 # 2 because of two sites
            
            #Etotal = Etot-EHF_local
            PhiDMFT = 2*TrSigmaG
            Etotal_DMFT = Etot - EHF_local + PhiDMFT  # 2 because of two sites
            #print 'E0=', Etot0, 'EH=', EH, 'EF=', EF, 'EHF=', 0.5*TrSG, 'EHF_local=', EHF_local, 'EDMFT=', PhiDMFT
            #print 'Etotal_DMFT=', Etotal_DMFT  # 2 because of two sites
            #
            #print 'N_local=', N_local
            #print 'Edc=', Edc
            #print 'EHF_local=', EHF_local
            #print 'GOOD=', sum([H0[i,i]*real(Nt2[i,i]) for i in range(len(H0))]) + 0.5*TrSG - U_local*N_local**2/4. * 2
            #print 'Tr(H0*G)+2/R=', sum([H0[i,i]*real(Nt2[i,i]) for i in range(len(H0))])+2/R
            #print 'E_{non-local-HF}-2/R=', 0.5*TrSG - U_local*N_local**2/4. * 2 - 2/R
            #print 'E_{non-local-HF}=', 0.5*TrSG-EHF_local
            #print 'Ediff=', abs(Ener[0]-p_Ener)
            #
            #print 'Nt2[In H2+ basis]='
            #for i in range(len(Evec)):
            #    for j in range(len(Evec)):
            #        print "%10.5f " % Nt2[i,j].real,
            #    print
            #print 'Nc[In H2+ basis]='
            #for i,(i1,i2) in enumerate(in12):
            #    if i2==0: print
            #    print "%10.5f "% Nc[i],
            #print
            
            print R, Etotal_DMFT, abs(Etotal_DMFT-p_Ener)
            savetxt(ctqmc.dir+'/'+'s_oo',[s_oo,[Edc]*len(s_oo)])
            
        savetxt('Nc.dat',Nc)
        
        if abs(Etotal_DMFT-p_Ener)<Etol: break
        p_Ener = Etotal_DMFT
