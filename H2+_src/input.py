from scipy import *


n_l_m = [(1,0,0),(2,1,0),(2,0,0)]

#nmax = 7
#n_l_m = [(1,0,0),(2,0,0)]
#for n in range(3,nmax+1):
#	n_l_m.append((n-1,1,0))
#	n_l_m.append((n,0,0))

nmax = 6
for n in range(3,nmax):
	for l in range(n):
		n_l_m.append((n,l,0))


print n_l_m
N_1ptl = len(n_l_m)
#print "nmax = ",nmax
print "# of (n,l,m) = ",N_1ptl
