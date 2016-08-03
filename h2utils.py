from scipy import *
import poly as polyx
from scipy import special

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
    maxk=len(aks)
    #Plm = polyx.Plmc(x,m,maxk)
    #mm = [aks[il]*Plm[il] for il in range(len(aks))]
    Plm = special.lpmn(m, maxk-1+m, x)[0][m]
    mm = [aks[il]*Plm[il+m] for il in range(len(aks))]
    return sum(mm)*exp(-q*x)

def Lambda(lam,m,p,Cks):
    "Function Lambda(lambda),part of the H2+ solution"
    alpha = m+0.0
    max_n = m+len(Cks)-1
    x = 2*p*(lam-1.)
    Lag = polyx.Laguerre_dirty(x,alpha,max_n)
    csum = sum([Cks[k]*Lag[k+m] for k in range(len(Cks))])  ### original
    #csum = sum([Cks[k]*special.genlaguerre(k,m)(x) for k in range(len(Cks))])  ### HERE Lag[k]->Lag[m+k]
    return sqrt((lam**2-1)**m)*exp(-x/2.)*csum

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

