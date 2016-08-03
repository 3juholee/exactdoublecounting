subroutine alp(l,m,p,A,R1, alp0,alp1,alp3)
  IMPLICIT NONE
  INTEGER, intent(in) :: l, m
  REAL*8, intent(in)  :: p, A, R1
  REAL*8, intent(out) :: alp0, alp1, alp3
  ! local variables
  REAL*8 :: q, q2
  INTEGER :: l2, l22, m2
  q = p
  l2 = l**2
  m2 = m**2
  l22 = 4*l2-1.
  q2 = q**2
  alp0 = (l2-m2)*(R1**2 - 4*q2*l2)/l22
  alp1 = A - q2 + l2-l
  alp3 = 8*p*l2 * (l2-m2)/l22
end subroutine alp
  
subroutine gam(k,m,p,A,R2, gam0, gam1, gam2, gam3)
  IMPLICIT NONE
  INTEGER, intent(in) :: k, m
  REAL*8, intent(in)  :: p, A, R2
  REAL*8, intent(out) :: gam0, gam1, gam2, gam3
  ! local variables
  REAL*8 :: p2, r2p
  p2 = p**2
  r2p = R2/(2*p)
  gam0 = (k+1)*(k+m+1)*(k+1-r2p)*(k+m+1-r2p)
  gam1 = A - p2 + R2 - (m+1)*(2*p+1-r2p) - 2*k*(k+m+2*p+1-r2p)
  gam2 = 2*p + (2*k+m+1)*(2+r2p/p)
  gam3 = (k+1)*(k+m+1)*(r2p/p)*(2*r2p-2*k-2-m)
end subroutine gam

subroutine BaberF(F0,p,A,m,Kmax,R1,start1)
  IMPLICIT NONE
  REAL*8, intent(out) :: F0
  REAL*8, intent(in)  :: p, A, R1, start1
  INTEGER, intent(in) :: m, Kmax
  !f2py integer optional, intent(in) :: start1=1e-10
  REAL*8  :: F2, F1, alp0, alp1, alp3
  INTEGER :: k
  F2 = 0
  F1 = start1
  do k=Kmax,0,-1
     CALL alp(m+k+1,m,p,A,R1, alp0,alp1,alp3)
     F0 = -alp1*F1 - alp0*F2
     F2 = F1
     F1 = F0
  end do
end subroutine BaberF
  
SUBROUTINE BaberD(F0,dFA0,dFp0,p,A,m,Kmax,R1,start1)
  IMPLICIT NONE
  REAL*8, intent(out) :: F0, dFA0, dFp0
  REAL*8, intent(in)  :: p, A, R1, start1
  INTEGER, intent(in) :: m, Kmax
  !f2py integer optional, intent(in) :: start1=1e-10
  REAL*8  :: F2, F1, alp0, alp1, alp3
  REAL*8  :: dFA1, dFA2, dFp1, dFp2
  INTEGER :: k
  F2 = 0
  F1 = start1
  dFA2 = 0
  dFA1 = 0
  dFp2 = 0
  dFp1 = 0
  do k=Kmax,0,-1
     CALL alp(m+k+1,m,p,A,R1, alp0,alp1,alp3)
     F0 = -alp1*F1 - alp0*F2
     dFA0 = -alp1*dFA1-alp0*dFA2-F1
     dFp0 = -alp1*dFp1-alp0*dFp2+2*p*F1+alp3*F2
     F2 = F1
     dFA2 = dFA1
     dFp2 = dFp1
     F1 = F0
     dFA1 = dFA0
     dFp1 = dFp0
  end do
end SUBROUTINE BaberD


SUBROUTINE HyllerF(G0,p,A,m,Kmax,R2,start1)
  IMPLICIT NONE
  REAL*8, intent(out) :: G0
  REAL*8, intent(in)  :: p, A, R2, start1
  INTEGER, intent(in) :: m, Kmax
  !f2py integer optional, intent(in) :: start1=1e-10
  REAL*8  :: G2, G1, gam0, gam1, gam2, gam3
  INTEGER :: k
  G2 = 0
  G1 = start1
  do k=Kmax,0,-1
     CALL gam(k,m,p,A,R2, gam0,gam1,gam2,gam3)     
     G0   = -gam1*G1   -gam0*G2
     G2 = G1
     G1 = G0
  enddo
end SUBROUTINE HyllerF

SUBROUTINE HyllerD(G0,dGA0,dGp0,p,A,m,Kmax,R2,start1)
  IMPLICIT NONE
  REAL*8, intent(out) :: G0, dGA0, dGp0
  REAL*8, intent(in)  :: p, A, R2, start1
  INTEGER, intent(in) :: m, Kmax
  !f2py integer optional, intent(in) :: start1=1e-10
  REAL*8  :: G2, G1, dGA2, dGA1, dGp2, dGp1, gam0, gam1, gam2, gam3
  INTEGER :: k
  G2 = 0
  G1 = start1
  dGA2 = 0
  dGA1 = 0
  dGp2 = 0
  dGp1 = 0
  do k=Kmax,0,-1
     CALL gam(k,m,p,A,R2, gam0,gam1,gam2,gam3)
     G0   = -gam1*G1   -gam0*G2
     dGA0 = -gam1*dGA1 -gam0*dGA2 - G1
     dGp0 = -gam1*dGp1 -gam0*dGp2 + gam2*G1 + gam3*G2
     G2 = G1
     dGA2 = dGA1
     dGp2 = dGp1
     G1 = G0
     dGA1 = dGA0
     dGp1 = dGp0
  end do
end SUBROUTINE HyllerD
