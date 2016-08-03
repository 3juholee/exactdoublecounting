from scipy import *
from scipy.misc import factorial2
from scipy import weave

def Plmc(x,m,maxk):
    Plm = zeros(maxk,dtype=float)
    Plm[0] = (-1)**m*factorial2(2*m-1)*sqrt((1-x**2)**m)
    Plm[1] = x*(2*m+1)*Plm[0]
    codep="""
      double xl = x;
      double Pm0 = Plm(0);
      double Pm1 = Plm(1);
      for (int l=m+1; l<m+maxk-1; l++){
          double c = l-m+1.;
          double beta = (2*l+1.)*xl/c;
          double gamma = (l+m)/c;
          double Pm2 = beta*Pm1-gamma*Pm0;
          Plm(l-m+1) = Pm2;
          Pm0 = Pm1;
          Pm1 = Pm2;
      }
    """
    weave.inline(codep, ['Plm', 'x', 'm', 'maxk'], type_converters=weave.converters.blitz, compiler = 'gcc')
    return Plm


def old_Laguerre_dirty(x,alpha,max_n):
    Lag = zeros(max_n+1, dtype=float)
    codel="""
       // Initialize the recursion process. //
       #include <iostream>
       using namespace std;
    
       double alpha_one_mx = (double) alpha + 1.0L - (double) x;
       if (max_n >=0 ){
       double ln2 = 1.0L;
       double ln1 = alpha_one_mx;
       Lag(0) = (double) ln2;
       if (max_n > 0){
       Lag(1) = (double) ln1;
       if (max_n > 1){
    
       // Calculate Ln(x) for n = 2,...,max_n //
    
       for (int k = 1; k < max_n; k++) {
          double beta = ( (double)(k + k) + alpha_one_mx );
          double gamma = (double) k + (double) alpha;
          double ln = (beta * ln1 - gamma * ln2) / (double)(k+1);
          Lag(k+1) = (double)ln;
          ln2 = ln1;
          ln1 = ln;
       }
       }
       }
       }
    
    """
    weave.inline(codel, ['Lag', 'x', 'alpha', 'max_n'], type_converters=weave.converters.blitz, compiler = 'gcc')
    return Lag
 
def Laguerre_dirty(x,alpha,max_n):
    Lag = zeros(max_n+1, dtype=float)
    codel="""
       // Initialize the recursion process. //
       #include <iostream>
       using namespace std;
    
       double alpha_one_mx = (double) alpha + 1.0L - (double) x;
       if (max_n >=0 ){
       double ln2 = 1.0L;
       double ln1 = alpha_one_mx;
       Lag(0) = (double) ln2;
       if (max_n > 0){
       Lag(1) = (double) ln1;
       if (max_n > 1){
    
       // Calculate Ln(x) for n = 2,...,max_n //
    
       for (int k = 1; k < max_n; k++) {
          double beta = ( (double)(k + k) + alpha_one_mx );
          double gamma = (double) k + (double) alpha;
          double ln = (beta * ln1 - gamma * ln2) / (double)(k+1);
          Lag(k+1) = (double)ln;
          ln2 = ln1;
          ln1 = ln;
       }
       }
       }
       }
    
    """
    weave.inline(codel, ['Lag', 'x', 'alpha', 'max_n'], type_converters=weave.converters.blitz, compiler = 'gcc')
    return Lag
    
   
if __name__ == '__main__':
    import sys
    from scipy import special
    
    max_n=100
    alpha=1.
    x=10.
    Lag = Laguerre_dirty(x,alpha,max_n)
    for n in range(max_n+1):
        print 'laguerre=', n, special.eval_genlaguerre(n,alpha,x), Lag[n]


    m=1
    max_l=30
    x=0.999
    
    Plm = Plmc(x,m,max_l)
    
    for il,l in enumerate(range(m,len(Plm)+m)):
        print 'Plm=', il,l,special.lpmv(m,l,x),Plm[il]
