!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program xtess
!
!    Test driver of "tess", the double precision Fortran subroutine to compute
!    the gravitational potential/acceleration vector/gradient tensor of a tesseroid.
!
!    Reference: 
!           Fukushima T (2018) Accurate computation of gravitational field of a tesseroid.
!                Journal of Geodesy 92(12):1371–1386, doi:10.1007/s00190-018-1126-2
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
PI=atan(1.d0)*4.d0
raddeg=PI/180.d0
delta=1.d-16
write (*,"(a15,1pe15.7)") "# delta=", delta
!
!    A sample tesseroid covering an eastern Himalaya: size 1 deg x 1 deg x 50 km
!      (Note) Mt. Everest: h = 8.848 km, phi = 27 deg 59 min, lat = 86 deg 55 min
!
R0=6380.d0
HB=-40.d0
HT=10.d0
PhiS=27.d0*raddeg
PhiN=28.d0*raddeg
FlamW=86.d0*raddeg
FlamE=87.d0*raddeg
!
!   Test evaluation points are along a radial straight line
!    passing through the geometrical center of the tesseroid
!
phi=27.5d0*raddeg
flam=86.5d0*raddeg
!
write (*,"(a15,1p3e15.7)") "# R0,HT,HB=",R0,HT,HB
write (*,"(a15,1p2e15.7)") "# PhiN,PhiS=",PhiN/raddeg,PhiS/raddeg
write (*,"(a15,1p2e15.7)") "# LamdaE,LamdaW=",FlamE/raddeg,FlamW/raddeg
write (*,"(a15,1p2e15.7)") "# phi,lambda=",phi/raddeg,flam/raddeg
!
write (*,"(a15,a25)") "# H","V"
write (*,"(a15,3a25)") "#","gPhi","gLambda","gH"
write (*,"(a15,3a25)") "#","GammaPhiPhi","GammaPhiLambda","GammaPhiH"
write (*,"(a15,3a25)") "#","GammaLambdaLambda","GammaLambdaH","GammaHH"
!
!    Compute the field along a radial line passing through the tesseroid
!
call cpu_time(time_begin)
dH=10.d0
do m=-10,29
    h=(dble(m)+0.5d0)*dH
!
    call tess(PhiN,PhiS,FlamE,FlamW,HT,HB,R0,delta, &
        phi,flam,h, &
        v,gPhi,gLam,gH,ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
    write (*,"(1pe15.7,1pe25.15)") h,v
    write (*,"(15x,1p3e25.15)") gPhi,gLam,gH
    write (*,"(15x,1p3e25.15)") ggPhiPhi,ggPhiLam,ggPhiH
    write (*,"(15x,1p3e25.15)") ggLamLam,ggLamH,ggHH
enddo
call cpu_time(time_end)
write(*,*) "The total time is",time_end-time_begin,"seconds."
!
stop
end program xtess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tess(PhiN0,PhiS0,FlamE0,FlamW0,HT0,HB0,R00,delta0,phi,flam,h, &
    v,gPhi,gLam,gH,ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
!    The double precision Fortran subroutine to compute
!    the gravitational potential, V, the gravitational acceleration vector, g,
!    and the gravity gradient tensor, Gamma, of a tesseroid.
!
!    Input parameters:
!      PhiN0: Latitude of the north end point of the tesseroid
!      PhiS0: Latitude of the south end point of the tesseroid, PhiS0 < PhiN0
!      FlamE0: Longitude of the east end point of the tesseroid
!      FlamW0: Longitude of the west end point of the tesseroid, FlamW0 < FlamE0
!      HT0: Height of the top end point of the tesseroid
!      HB0: Height of the bottom end point of the tesseroid, HB0 < HT0
!      R00: Radius of the reference spherical surface from which H is counted
!      delta0: the relative error tolerance (typical value is 1d-9 or 1d-12)
!    Input variables:
!      phi: Latitude of the evaluation point
!      flam: Longitude of the evaluation point
!      h: Height from the reference sphere of the evaluation point
!        (Note) We adopted "km" as the unit of h as well as HT0 and HB0.
!               But any unit is OK if R00 is expressed in the same unit.
!    Output variables:
!      v: Normalized gravitational potential at the evaluation point
!        (Note) normalization constant: V0 = G rho R00^2
!      gPhi: Latitude component of normalized gravitational acceleration vector
!      gLam: Longitude component of normalized gravitational acceleration vector
!      gH: Height component of normalized gravitational acceleration vector
!        (Note) normalization constant: g0 = G rho R00
!      ggPhiPhi: Latitude-latitude component of normalized Gamma
!      ggPhiLam: Latitude-longitude component of normalized Gamma
!      ggPhiH: Latitude-height component of normalized Gamma
!      ggLamLam: Longitude-longitude component of normalized Gamma
!      ggLamH: Longitude-height component of normalized Gamma
!      ggHH: Height-height component of normalized Gamma
!        (Note) normalization constant: Gamma0 = G rho
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "vtess" as an external function name
!
external vtess
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
!
!    Store the input variables into common blocks
!    so that the subprograms can refer them indirectly
!
PhiN=PhiN0
PhiS=PhiS0
FlamE=FlamE0
FlamW=FlamW0
HT=HT0
HB=HB0
R0=R00
!
!    Calculate parameters of numerical integration and differentiation
!
eps=delta0
delta=delta0
delta1=delta**(1.d0/3.d0)
delta2=sqrt(sqrt(delta))
!
!    Compute tesseroid-dependent constants beforehand
!
deltaH=HT-HB
deltaPhi=PhiN-PhiS
deltaFlam=FlamE-FlamW
beta=deltaH/R0
RC=R0+(HT+HB)*0.5d0
PhiC=(PhiN+PhiS)*0.5d0
FlamC=(FlamE+FlamW)*0.5d0
sinPhiC=sin(PhiC)
cosPhiC=cos(PhiC)
!
!    Compute evaluation-point-dependent constants beforehand
!
r=R0+h
sinphi=sin(phi)
cosphi=cos(phi)
cosdlam=cos(flam-FlamC)
scale0=0.5d0*sqrt(deltaH*deltaH+RC*RC*(deltaPhi*deltaPhi &
    +cosPhiC*cosPhiC*deltaFlam*deltaFlam))
scale1=sqrt(r*r+RC*RC-2.d0*RC*r*(sinphi*sinPhiC+cosphi*cosPhiC*cosdlam))
scale=max(scale0,scale1)
d1h=scale*delta1
! add d1phi, d1lam
d1phi=d1h/r
d1lam=d1h/(r*cosphi)
d2h=scale*delta2
! add d2phi, d2lam
d2phi=d2h/r
d2lam=d2h/(r*cosphi)
!
!    Compute the gravitational potential by numerical integration
!
v=vtess(phi,flam,h)
!
!    Compute the gravitational acceleration vector by numerical differentiation
!
call gtess(vtess,phi,flam,h,v,gPhi,gLam,gH)
!
!    Compute the gravity gradient tensor by numerical differentiation
!
call ggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function vtess(phi,flam,h)
!
!    The double precision Fortran function to compute
!    the gravitational potential by the conditional split quadrature method
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "stess" as an external function name
!
external stess
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /argV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Store the input variables into common blocks
!    so that the subprograms can refer them indirectly
!
phi0=phi
flam0=flam
h0=h
!
!    Compute tesseroid-dependent constants beforehand
!
xiN=PhiN-phi
xiS=PhiS-phi
etaE=FlamE-flam
etaW=FlamW-flam
!
!    Compute evaluation-point-dependent constants beforehand
!
alpha=1.d0+h/R0
gamma=cos(phi)
zetaB=(HB-h)/R0
zetaT=zetaB+beta
!
!    Conditional split quadrature of latitude line integration 
!
if(xiS.lt.0.d0.and.xiN.gt.0.d0) then
    call dqde(stess,xiS,0.d0,eps,v1)
    call dqde(stess,0.d0,xiN,eps,v2)
    v=v1+v2
else
    call dqde(stess,xiS,xiN,eps,v)
endif
vtess=v
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function stess(xi)
!
!    The double precision Fortran function to compute
!    the latitude line integral requried in "vtess"
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "vkern" as an external function name
!
external vkern
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /argV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Compute latitude-dependent constants beforehand
!
sigma=cos(phi0+xi)
fmu=sin(xi*0.5d0)
!
!    Conditional split quadrature of longitude line integration 
!
if(etaW.lt.0.d0.and.etaE.gt.0.d0) then
    call dqde1(vkern,etaW,0.d0,eps,s1)
    call dqde1(vkern,0.d0,etaE,eps,s2)
    s=s1+s2
else
    call dqde1(vkern,etaW,etaE,eps,s)
endif
stess=s
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function vkern(eta)
!
!    The double precision Fortran function to compute
!    the kernel function expressed in a cancellation-free form
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Cancellation-free computation of the kernel function
!
fnu=sin(eta*0.5d0)
B=2.d0*alpha*(fmu*fmu+gamma*sigma*fnu*fnu)
A=B*(2.d0*alpha-B)
C=(alpha-B)*(alpha-B)-A*0.5d0
D=0.5d0*(zetaT+4.d0*alpha-3.d0*B)
ST=sqrt(A+(B+zetaT)*(B+zetaT))
SB=sqrt(A+(B+zetaB)*(B+zetaB))
T=(zetaT+zetaB+2.d0*B)/(ST+SB)
if(zetaB+B.lt.0.d0) then
    fl=dlog1p((1.d0+T)*beta*(-zetaB-B+SB)/A)
else
    fl=dlog1p((1.d0+T)*beta/(zetaB+B+SB))
endif
vkern=sigma*(C*fl+beta*(D*T+SB*0.5d0))
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gtess(vtess,phi,flam,h,v,gPhi,gLam,gH)
!
!    Compute the gravitational acceleration vector by the conditional switch
!    of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
!
!    Compute conmponent-independent variables beorehand
!
r=R0+h
c=cos(phi)
!
!   Compute gPhi
!
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d1phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d1phi.lt.PhiN)) then
        gPhi=(-vtess(phi+2.d0*d1phi,flam,h)+4.d0*vtess(phi+d1phi,flam,h)-3.d0*v)/(2.d0*d1phi)
    elseif((phi.lt.PhiS.and.phi+d1phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d1phi.gt.PhiN)) then
        gPhi=(vtess(phi-2.d0*d1phi,flam,h)-4.d0*vtess(phi-d1phi,flam,h)+3.d0*v)/(2.d0*d1phi)
    else
        gPhi=(vtess(phi+d1phi,flam,h)-vtess(phi-d1phi,flam,h))/(2.d0*d1phi)
    endif
else
    gPhi=(vtess(phi+d1phi,flam,h)-vtess(phi-d1phi,flam,h))/(2.d0*d1phi)
endif
gPhi=gPhi/r
!
!   Compute gLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d1lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d1lam.lt.FlamE)) then
        gLam=(-vtess(phi,flam+2.d0*d1lam,h)+4.d0*vtess(phi,flam+d1lam,h)-3.d0*v)/(2.d0*d1lam)
    elseif((flam.lt.FlamW.and.flam+d1lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d1lam.gt.FlamE)) then
        gLam=(vtess(phi,flam-2.d0*d1lam,h)-4.d0*vtess(phi,flam-d1lam,h)+3.d0*v)/(2.d0*d1lam)
    else
        gLam=(vtess(phi,flam+d1lam,h)-vtess(phi,flam-d1lam,h))/(2.d0*d1lam)
    endif
else
    gLam=(vtess(phi,flam+d1lam,h)-vtess(phi,flam-d1lam,h))/(2.d0*d1lam)
endif
gLam=gLam/(r*c)
!
!   Compute gH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d1h.lt.HB).or.(h.gt.HT.and.h-d1h.lt.HT)) then
        gH=(-vtess(phi,flam,h+2.d0*d1h)+4.d0*vtess(phi,flam,h+d1h)-3.d0*v)/(2.d0*d1h)
    elseif((h.lt.HB.and.h+d1h.gt.HB).or.(h.lt.HT.and.h+d1h.gt.HT)) then
        gH=(vtess(phi,flam,h-2.d0*d1h)-4.d0*vtess(phi,flam,h-d1h)+3.d0*v)/(2.d0*d1h)
    else
        gH=(vtess(phi,flam,h+d1h)-vtess(phi,flam,h-d1h))/(2.d0*d1h)
    endif
else
    gH=(vtess(phi,flam,h+d1h)-vtess(phi,flam,h-d1h))/(2.d0*d1h)
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
!    Compute the gravity gradient tensor by the conditional switch
!    of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h
!
!    Compute conmponent-independent variables beorehand
!
r=R0+h
c=cos(phi)
t=tan(phi)
!
!   Compute ggPhiPhi
!
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)) then
        ggPhiPhi=(-vtess(phi+3.d0*d2phi,flam,h)+4.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -5.d0*vtess(phi+d2phi,flam,h)+2.d0*v)/(d2phi*d2phi)
    elseif((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)) then
        ggPhiPhi=(-vtess(phi-3.d0*d2phi,flam,h)+4.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -5.d0*vtess(phi-d2phi,flam,h)+2.d0*v)/(d2phi*d2phi)
    else
        ggPhiPhi=(vtess(phi+d2phi,flam,h)-2.d0*v+vtess(phi-d2phi,flam,h))/(d2phi*d2phi)
    endif
else
    ggPhiPhi=(vtess(phi+d2phi,flam,h)-2.d0*v+vtess(phi-d2phi,flam,h))/(d2phi*d2phi)
endif
ggPhiPhi=ggPhiPhi/(r*r)+gH/r
!
!   Compute ggLamLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)) then
        ggLamLam=(-vtess(phi,flam+3.d0*d2lam,h)+4.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -5.d0*vtess(phi,flam+d2lam,h)+2.d0*v)/(d2lam*d2lam)
    elseif((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)) then
        ggLamLam=(-vtess(phi,flam-3.d0*d2lam,h)+4.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -5.d0*vtess(phi,flam-d2lam,h)+2.d0*v)/(d2lam*d2lam)
    else
        ggLamLam=(vtess(phi,flam+d2lam,h)-2.d0*v+vtess(phi,flam-d2lam,h))/(d2lam*d2lam)
    endif
else
    ggLamLam=(vtess(phi,flam+d2lam,h)-2.d0*v+vtess(phi,flam-d2lam,h))/(d2lam*d2lam)
endif
ggLamLam=ggLamLam/(r*r*c*c)+(gH-gPhi*t)/r
!
!   Compute ggHH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)) then
        ggHH=(-vtess(phi,flam,h+3.d0*d2h)+4.d0*vtess(phi,flam,h+2.d0*d2h) &
            -5.d0*vtess(phi,flam,h+d2h)+2.d0*v)/(d2h*d2h)
    elseif((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)) then
        ggHH=(-vtess(phi,flam,h-3.d0*d2h)+4.d0*vtess(phi,flam,h-2.d0*d2h) &
            -5.d0*vtess(phi,flam,h-d2h)+2.d0*v)/(d2h*d2h)
    else
        ggHH=(vtess(phi,flam,h+d2h)-2.d0*v+vtess(phi,flam,h-d2h))/(d2h*d2h)
    endif
else
    ggHH=(vtess(phi,flam,h+d2h)-2.d0*v+vtess(phi,flam,h-d2h))/(d2h*d2h)
endif
!
!   Compute ggPhiLam
!
if(h.ge.HB.and.h.le.HT) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=(vtess(phi+2.d0*d2phi,flam+2.d0*d2lam,h) &
            -4.d0*vtess(phi+2.d0*d2phi,flam+d2lam,h)+3.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi+d2phi,flam+2.d0*d2lam,h)+16.d0*vtess(phi+d2phi,flam+d2lam,h) &
            -12.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam+d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=-(vtess(phi+2.d0*d2phi,flam-2.d0*d2lam,h) &
            -4.d0*vtess(phi+2.d0*d2phi,flam-d2lam,h)+3.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi+d2phi,flam-2.d0*d2lam,h)+16.d0*vtess(phi+d2phi,flam-d2lam,h) &
            -12.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam-d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=-(vtess(phi-2.d0*d2phi,flam+2.d0*d2lam,h) &
            -4.d0*vtess(phi-2.d0*d2phi,flam+d2lam,h)+3.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi-d2phi,flam+2.d0*d2lam,h)+16.d0*vtess(phi-d2phi,flam+d2lam,h) &
            -12.d0*vtess(phi-d2phi,flam,h)+3.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam+d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=(vtess(phi-2.d0*d2phi,flam-2.d0*d2lam,h) &
            -4.d0*vtess(phi-2.d0*d2phi,flam-d2lam,h)+3.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi-d2phi,flam-2.d0*d2lam,h)+16.d0*vtess(phi-d2phi,flam-d2lam,h) &
            -12.d0*vtess(phi-d2phi,flam,h)+3.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam-d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(-vtess(phi+2.d0*d2phi,flam+d2lam,h)+vtess(phi+2.d0*d2phi,flam-d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam+d2lam,h)-4.d0*vtess(phi+d2phi,flam-d2lam,h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(vtess(phi-2.d0*d2phi,flam+d2lam,h)-vtess(phi-2.d0*d2phi,flam-d2lam,h) &
            -4.d0*vtess(phi-d2phi,flam+d2lam,h)+4.d0*vtess(phi-d2phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam+d2lam,h)-3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or. &
        (flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=(-vtess(phi+d2phi,flam+2.d0*d2lam,h)+vtess(phi-d2phi,flam+2.d0*d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam+d2lam,h)-4.d0*vtess(phi-d2phi,flam+d2lam,h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2lam)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or. &
        (flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=-(-vtess(phi+d2phi,flam-2.d0*d2lam,h)+vtess(phi-d2phi,flam-2.d0*d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam-d2lam,h)-4.d0*vtess(phi-d2phi,flam-d2lam,h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2lam)
    else
        ggPhiLam=(vtess(phi+d2phi,flam+d2lam,h)-vtess(phi-d2phi,flam+d2lam,h) &
            -vtess(phi+d2phi,flam-d2lam,h)+vtess(phi-d2phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    endif
else
    ggPhiLam=(vtess(phi+d2phi,flam+d2lam,h)-vtess(phi-d2phi,flam+d2lam,h) &
        -vtess(phi+d2phi,flam-d2lam,h)+vtess(phi-d2phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
endif
ggPhiLam=ggPhiLam/(r*r*c)+gLam*t/r
!
!   Compute ggPhiH
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=(vtess(phi+2.d0*d2phi,flam,h+2.d0*d2h)-4.d0*vtess(phi+2.d0*d2phi,flam,h+d2h) &
            +3.d0*vtess(phi+2.d0*d2phi,flam,h)-4.d0*vtess(phi+d2phi,flam,h+2.d0*d2h) &
            +16.d0*vtess(phi+d2phi,flam,h+d2h)-12.d0*vtess(phi+d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=-(vtess(phi+2.d0*d2phi,flam,h-2.d0*d2h)-4.d0*vtess(phi+2.d0*d2phi,flam,h-d2h) &
            +3.d0*vtess(phi+2.d0*d2phi,flam,h)-4.d0*vtess(phi+d2phi,flam,h-2.d0*d2h) &
            +16.d0*vtess(phi+d2phi,flam,h-d2h)-12.d0*vtess(phi+d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=-(vtess(phi-2.d0*d2phi,flam,h+2.d0*d2h)-4.d0*vtess(phi-2.d0*d2phi,flam,h+d2h) &
            +3.d0*vtess(phi-2.d0*d2phi,flam,h)-4.d0*vtess(phi-d2phi,flam,h+2.d0*d2h) &
            +16.d0*vtess(phi-d2phi,flam,h+d2h)-12.d0*vtess(phi-d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=(vtess(phi-2.d0*d2phi,flam,h-2.d0*d2h)-4.d0*vtess(phi-2.d0*d2phi,flam,h-d2h) &
            +3.d0*vtess(phi-2.d0*d2phi,flam,h)-4.d0*vtess(phi-d2phi,flam,h-2.d0*d2h) &
            +16.d0*vtess(phi-d2phi,flam,h-d2h)-12.d0*vtess(phi-d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=(-vtess(phi+2.d0*d2phi,flam,h+d2h)+vtess(phi+2.d0*d2phi,flam,h-d2h) &
            +4.d0*vtess(phi+d2phi,flam,h+d2h)-4.d0*vtess(phi+d2phi,flam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=-(-vtess(phi-2.d0*d2phi,flam,h+d2h)+vtess(phi-2.d0*d2phi,flam,h-d2h) &
            +4.d0*vtess(phi-d2phi,flam,h+d2h)-4.d0*vtess(phi-d2phi,flam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=(-vtess(phi+d2phi,flam,h+2.d0*d2h)+vtess(phi-d2phi,flam,h+2.d0*d2h) &
            +4.d0*vtess(phi+d2phi,flam,h+d2h)-4.d0*vtess(phi-d2phi,flam,h+d2h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=-(-vtess(phi+d2phi,flam,h-2.d0*d2h)+vtess(phi-d2phi,flam,h-2.d0*d2h) &
            +4.d0*vtess(phi+d2phi,flam,h-d2h)-4.d0*vtess(phi-d2phi,flam,h-d2h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2h)
    else
        ggPhiH=(vtess(phi+d2phi,flam,h+d2h)-vtess(phi-d2phi,flam,h+d2h) &
            -vtess(phi+d2phi,flam,h-d2h)+vtess(phi-d2phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    endif
else
    ggPhiH=(vtess(phi+d2phi,flam,h+d2h)-vtess(phi-d2phi,flam,h+d2h) &
        -vtess(phi+d2phi,flam,h-d2h)+vtess(phi-d2phi,flam,h-d2h))/(4.d0*d2phi*d2h)
endif
ggPhiH=ggPhiH/r-gPhi/r
!
!   Compute ggLamH
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
    if(((flam.gt.FlamW.and.flam-d2flam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2flam.lt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=(vtess(phi,flam+2.d0*d2lam,h+2.d0*d2h)-4.d0*vtess(phi,flam+2.d0*d2lam,h+d2h) &
            +3.d0*vtess(phi,flam+2.d0*d2lam,h)-4.d0*vtess(phi,flam+d2lam,h+2.d0*d2h) &
            +16.d0*vtess(phi,flam+d2lam,h+d2h)-12.d0*vtess(phi,flam+d2lam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=-(vtess(phi,flam+2.d0*d2lam,h-2.d0*d2h)-4.d0*vtess(phi,flam+2.d0*d2lam,h-d2h) &
            +3.d0*vtess(phi,flam+2.d0*d2lam,h)-4.d0*vtess(phi,flam+d2lam,h-2.d0*d2h) &
            +16.d0*vtess(phi,flam+d2lam,h-d2h)-12.d0*vtess(phi,flam+d2lam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=-(vtess(phi,flam-2.d0*d2lam,h+2.d0*d2h)-4.d0*vtess(phi,flam-2.d0*d2lam,h+d2h) &
            +3.d0*vtess(phi,flam-2.d0*d2lam,h)-4.d0*vtess(phi,flam-d2lam,h+2.d0*d2h) &
            +16.d0*vtess(phi,flam-d2lam,h+d2h)-12.d0*vtess(phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=(vtess(phi,flam-2.d0*d2lam,h-2.d0*d2h)-4.d0*vtess(phi,flam-2.d0*d2lam,h-d2h) &
            +3.d0*vtess(phi,flam-2.d0*d2lam,h)-4.d0*vtess(phi,flam-d2lam,h-2.d0*d2h) &
            +16.d0*vtess(phi,flam-d2lam,h-d2h)-12.d0*vtess(phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=(-vtess(phi,flam+2.d0*d2lam,h+d2h)+vtess(phi,flam+2.d0*d2lam,h-d2h) &
            +4.d0*vtess(phi,flam+d2lam,h+d2h)-4.d0*vtess(phi,flam+d2lam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=-(-vtess(phi,flam-2.d0*d2lam,h+d2h)+vtess(phi,flam-2.d0*d2lam,h-d2h) &
            +4.d0*vtess(phi,flam-d2lam,h+d2h)-4.d0*vtess(phi,flam-d2lam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2lam*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=(-vtess(phi,flam+d2lam,h+2.d0*d2h)+vtess(phi,flam-d2lam,h+2.d0*d2h) &
            +4.d0*vtess(phi,flam+d2lam,h+d2h)-4.d0*vtess(phi,flam-d2lam,h+d2h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2lam*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=-(-vtess(phi,flam+d2lam,h-2.d0*d2h)+vtess(phi,flam-d2lam,h-2.d0*d2h) &
            +4.d0*vtess(phi,flam+d2lam,h-d2h)-4.d0*vtess(phi,flam-d2lam,h-d2h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2lam*d2h)
    else
        ggLamH=(vtess(phi,flam+d2lam,h+d2h)-vtess(phi,flam-d2lam,h+d2h) &
            -vtess(phi,flam+d2lam,h-d2h)+vtess(phi,flam-d2lam,h-d2h))/(4.d0*d2lam*d2h)
    endif
else
    ggLamH=(vtess(phi,flam+d2lam,h+d2h)-vtess(phi,flam-d2lam,h+d2h) &
        -vtess(phi,flam+d2lam,h-d2h)+vtess(phi,flam-d2lam,h-d2h))/(4.d0*d2lam*d2h)
endif
ggLamH=ggLamH/(r*c)-gLam/r
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dqde(f,a,b,delta0,s)
!
!   The subroutine to compute a given integral by the double precision
!    double exponential (DE) quadrature rule.
!
!    Reference: Fukushima (2017) Precise and fast computation of gravitational
!      field of general finite body and its application to gravitational study
!      of asteroid Eros. Astron J. in printing
!
!    Input function/variables:
!      f(t): a function to be integrated 
!      a: the lower end point of the integration interval
!      b: the upper end point of the integration interval
!      delta0: the relative error tolerance (typically 1d-15)
!    Output variable:
!      s: the definite integral, s=int_a^b f(t) dt
!        (Note) "errd" is the latest relative error estimate o "s"
!
real*8 f,a,b,delta0,s
integer MMAX
real*8 SAFETY,TWOPI
parameter (MMAX=1024,SAFETY=10.d0)
parameter (TWOPI=6.2831853071795865d0)
integer m
real*8 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*8 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*8 errh,errd
delta=max(1.d-16,min(0.01d0,delta0))
if(delta.gt.1.d-4) then
    hmax=2.75d0+(-log10(delta))*0.75d0
elseif(delta.gt.1.d-6) then
    hmax=3.75d0+(-log10(delta))*0.5d0
elseif(delta.gt.1.d-10) then
    hmax=5.25d0+(-log10(delta))*0.25d0
else
    hmax=6.5d0+(-log10(delta))*0.125d0
endif
deltax=SAFETY*delta
factor=1.d0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.d0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5d0)*(bma*0.25d0)
s=sr*(2.d0*TWOPI)
err=abs(s)*deltat
h=2.d0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5d0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.d0/(1.d0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.d0-xw)
            fa=f(a+xa)*wg
            fb=f(b-xa)*wg
            fapfb=fa+fb
            sr=sr+fapfb
            s=s+fapfb*eppem
            errt=(abs(fa)+abs(fb))*eppem
            if(m.eq.1) err=err+errt*deltat
            ep=ep*eph
            em=em*emh
            if(errt.gt.err.or.xw.gt.deltah) goto 3
        t=t+h
        if(t.lt.h0) goto 2
    if(m.eq.1) then
        errh=(err/deltat)*deltah*h0
        errd=1.d0+2.d0*errh
    else
        errd=h*(abs(s-2.d0*sprev)+4.d0*abs(sr-2.d0*srprev))
    endif
    h=h*0.5d0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(dqde) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dqde1(f,a,b,delta0,s)
!
!    The first copy of "dqde" to be used in computing the internal
!    line integral called by "dqde"
!
real*8 f,a,b,delta0,s
integer MMAX
real*8 SAFETY,TWOPI
parameter (MMAX=1024,SAFETY=10.d0)
parameter (TWOPI=6.2831853071795865d0)
integer m
real*8 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*8 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*8 errh,errd
delta=max(1.d-16,min(0.01d0,delta0))
if(delta.gt.1.d-4) then
    hmax=2.75d0+(-log10(delta))*0.75d0
elseif(delta.gt.1.d-6) then
    hmax=3.75d0+(-log10(delta))*0.5d0
elseif(delta.gt.1.d-10) then
    hmax=5.25d0+(-log10(delta))*0.25d0
else
    hmax=6.5d0+(-log10(delta))*0.125d0
endif
deltax=SAFETY*delta
factor=1.d0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.d0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5d0)*(bma*0.25d0)
s=sr*(2.d0*TWOPI)
err=abs(s)*deltat
h=2.d0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5d0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.d0/(1.d0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.d0-xw)
            fa=f(a+xa)*wg
            fb=f(b-xa)*wg
            fapfb=fa+fb
            sr=sr+fapfb
            s=s+fapfb*eppem
            errt=(abs(fa)+abs(fb))*eppem
            if(m.eq.1) err=err+errt*deltat
            ep=ep*eph
            em=em*emh
            if(errt.gt.err.or.xw.gt.deltah) goto 3
        t=t+h
        if(t.lt.h0) goto 2
    if(m.eq.1) then
        errh=(err/deltat)*deltah*h0
        errd=1.d0+2.d0*errh
    else
        errd=h*(abs(s-2.d0*sprev)+4.d0*abs(sr-2.d0*srprev))
    endif
    h=h*0.5d0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(dqde1) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function dlog1p(x0)
!
!    The double precision Fortran function to compute precisely log1p(x),
!    a special logarithm function defined as log1p(x) = ln(1+x),
!    by the (7,6) rational minimax approximation evaluated by Estrin's scheme
!
!    Reference: Fukushima (2017) Precise and fast computation of gravitational
!      field of general finite body and its application to gravitational study
!      of asteroid Eros. Astron J. in printing
!
real*8 x0,x1,x2,x4,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6
parameter (a1=0.99999999999999999405d0,a2=2.45235728562912886048d0)
parameter (a3=2.17053627298972253249d0,a4=0.83928994566440838378d0)
parameter (a5=0.13520496594993836479d0,a6=0.00682631751459270270d0)
parameter (a7=0.00002291289324181940d0)
parameter (b1=2.95235728562912599232d0,b2=3.31338158247117791600d0)
parameter (b3=1.76186164168333482938d0,b4=0.44976458082070468584d0)
parameter (b5=0.04896199808811261680d0,b6=0.00157389087429218809d0)
if(x0.lt.-0.5d0.or.x0.gt.1.d0) then
    dlog1p=log(1.d0+x0)
elseif(x0.lt.0.d0) then
    x1=-x0/(1.d0+x0)
    x2=x1*x1; x4=x2*x2
    dlog1p=-x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
        /(((1.d0+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
else
    x1=x0
    x2=x1*x1; x4=x2*x2
    dlog1p=x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
        /(((1.d0+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Output of "xtess"
!      (Note) It took 3.568 s at a PC with an Intel Core i5-10400 running at 2.90 GHz
!      The ifort version is
!           Intel(R) Fortran Intel(R) 64 Compiler Classic for applications running on 
!            Intel(R) 64, Version 2021.4.0 Build 20210910_000000
!            Copyright (C) 1985-2021 Intel Corporation.  All rights reserved.
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        # delta=  1.0000000E-16
!     # R0,HT,HB=  6.3800000E+03  1.0000000E+01 -4.0000000E+01
!    # PhiN,PhiS=  2.8000000E+01  2.7000000E+01
! # LamdaE,LamdaW  8.7000000E+01  8.6000000E+01
!   # phi,lambda=  2.7500000E+01  8.6500000E+01
!             # H                        V
!               #                     gPhi                  gLambda                       gH
!               #              GammaPhiPhi           GammaPhiLambda                GammaPhiH
!               #        GammaLambdaLambda             GammaLambdaH                  GammaHH
!  -9.5000000E+01    1.522385430859932E-04
!                   -6.913732182770679E-10    0.000000000000000E+00    1.564130608898578E-06
!                   -1.387341662716454E-08    0.000000000000000E+00   -1.613894697029786E-11
!                   -1.505917753693464E-08    0.000000000000000E+00    2.893259481144160E-08
!  -8.5000000E+01    1.694526871210555E-04
!                   -8.747066001775616E-10    0.000000000000000E+00    1.892198142307468E-06
!                   -1.760580950749764E-08    0.000000000000000E+00   -2.066058851304799E-11
!                   -1.943435020150582E-08    0.000000000000000E+00    3.704016020364891E-08
!  -7.5000000E+01    1.903905934236350E-04
!                   -1.106672542712507E-09    0.000000000000000E+00    2.312822513933105E-06
!                   -2.233996420953945E-08    0.000000000000000E+00   -2.578804971985919E-11
!                   -2.516840640266663E-08    0.000000000000000E+00    4.750836836462066E-08
!  -6.5000000E+01    2.161011857153642E-04
!                   -1.389893134484920E-09    0.000000000000000E+00    2.851074036457655E-06
!                   -2.812502336397519E-08    0.000000000000000E+00   -3.071408279579221E-11
!                   -3.244537518351365E-08    0.000000000000000E+00    6.057039583336358E-08
!  -5.5000000E+01    2.478880058505530E-04
!                   -1.715013458072009E-09    0.000000000000000E+00    3.531846611828729E-06
!                   -3.475500622664502E-08    0.000000000000000E+00   -3.383415931884386E-11
!                   -4.112692458673271E-08    0.000000000000000E+00    7.588192820201063E-08
!  -4.5000000E+01    2.872708844104648E-04
!                   -2.053259639166731E-09    0.000000000000000E+00    4.371789329760953E-06
!                   -4.161825034453409E-08    0.000000000000000E+00   -3.297612891861420E-11
!                   -5.045690145966183E-08    0.000000000000000E+00    9.207514074757049E-08
!  -3.5000000E+01    3.319910086311147E-04
!                   -2.355567528758799E-09    0.000000000000000E+00    3.825615688241572E-06
!                   -4.793266758434474E-08    0.000000000000000E+00   -2.649421427519882E-11
!                   -5.920230924283715E-08    0.000000000000000E+00   -2.015876831867747E-07
!  -2.5000000E+01    3.603747131612369E-04
!                   -2.564755213201883E-09    0.000000000000000E+00    1.869421744025771E-06
!                   -5.251553452557477E-08    0.000000000000000E+00   -1.457160314363856E-11
!                   -6.555716477278904E-08    0.000000000000000E+00   -1.906499670213037E-07
!  -1.5000000E+01    3.696350462514114E-04
!                   -2.636222132600265E-09    0.000000000000000E+00   -1.064260798473823E-08
!                   -5.421371523781516E-08    0.000000000000000E+00    5.471010114351434E-13
!                   -6.787327496614387E-08    0.000000000000000E+00   -1.866356603343685E-07
!  -5.0000000E+00    3.601669342163731E-04
!                   -2.554314977920890E-09    0.000000000000000E+00   -1.889139936098946E-06
!                   -5.272632315536577E-08    0.000000000000000E+00    1.551673604698940E-11
!                   -6.568213752785291E-08    0.000000000000000E+00   -1.903142199554053E-07
!   5.0000000E+00    3.316111352803727E-04
!                   -2.337440063229249E-09    0.000000000000000E+00   -3.839348173993404E-06
!                   -4.846251165299864E-08    0.000000000000000E+00    2.705997748246749E-11
!                   -5.959167367262047E-08    0.000000000000000E+00   -2.006684780520396E-07
!   1.5000000E+01    2.868101440205006E-04
!                   -2.031662182754119E-09    0.000000000000000E+00   -4.374218164388472E-06
!                   -4.215634363810153E-08    0.000000000000000E+00    3.311074245470626E-11
!                   -5.084513496897725E-08    0.000000000000000E+00    9.300147090783754E-08
!   2.5000000E+01    2.474382990881564E-04
!                   -1.693806956102326E-09    0.000000000000000E+00   -3.528159772641221E-06
!                   -3.499501983026314E-08    0.000000000000000E+00    3.364904182370082E-11
!                   -4.124207840873123E-08    0.000000000000000E+00    7.623710370309283E-08
!   3.5000000E+01    2.157004300813909E-04
!                   -1.371418617690112E-09    0.000000000000000E+00   -2.845439875551469E-06
!                   -2.820499777133556E-08    0.000000000000000E+00    3.038054237659985E-11
!                   -3.243694585604362E-08    0.000000000000000E+00    6.064194108788116E-08
!   4.5000000E+01    1.900473207915234E-04
!                   -1.091721125483427E-09    0.000000000000000E+00   -2.307151583258427E-06
!                   -2.234592357243641E-08    0.000000000000000E+00    2.543375054579544E-11
!                   -2.511765420452311E-08    0.000000000000000E+00    4.746358413013454E-08
!   5.5000000E+01    1.691630798625039E-04
!                   -8.631048686893352E-10    0.000000000000000E+00   -1.887193033025766E-06
!                   -1.758321987787706E-08    0.000000000000000E+00    2.035095600256552E-11
!                   -1.937716647911629E-08    0.000000000000000E+00    3.696038720953256E-08
!   6.5000000E+01    1.519948957924415E-04
!                   -6.825568352191813E-10    0.000000000000000E+00   -1.559945026069056E-06
!                   -1.384333888277125E-08    0.000000000000000E+00    1.589237553727608E-11
!                   -1.500835554828129E-08    0.000000000000000E+00    2.885169681223171E-08
!   7.5000000E+01    1.377280082166444E-04
!                   -5.424089048682221E-10    0.000000000000000E+00   -1.303706863264856E-06
!                   -1.095668253283118E-08    0.000000000000000E+00    1.227887763567634E-11
!                   -1.172281702557548E-08    0.000000000000000E+00    2.267949631251105E-08
!   8.5000000E+01    1.257414395776617E-04
!                   -4.342889686494518E-10    0.000000000000000E+00   -1.101402778726971E-06
!                   -8.741568381379833E-09    0.000000000000000E+00    9.468292438805019E-12
!                   -9.254182417357096E-09    0.000000000000000E+00    1.799575571025114E-08
!   9.5000000E+01    1.155637184360951E-04
!                   -3.508161513851574E-10    0.000000000000000E+00   -9.400780632304917E-07
!                   -7.039686218393687E-09    0.000000000000000E+00    7.324354890423470E-12
!                   -7.389355216278842E-09    0.000000000000000E+00    1.442904793831886E-08
!   1.0500000E+02    1.068357959386216E-04
!                   -2.860576710892528E-10    0.000000000000000E+00   -8.100601082309440E-07
!                   -5.725034047993713E-09    0.000000000000000E+00    5.702126053128496E-12
!                   -5.968287512329931E-09    0.000000000000000E+00    1.169332748275438E-08
!   1.1500000E+02    9.928233337092832E-05
!                   -2.354546021701575E-10    0.000000000000000E+00   -7.041574277348270E-07
!                   -4.701624500135331E-09    0.000000000000000E+00    4.474584025759874E-12
!                   -4.874110373984178E-09    0.000000000000000E+00    9.575735816747436E-09
!   1.2500000E+02    9.269027708860816E-05
!                   -1.955797189911766E-10    0.000000000000000E+00   -6.170099768474632E-07
!                   -3.897812949189702E-09    0.000000000000000E+00    3.542400537102165E-12
!                   -4.022352420927403E-09    0.000000000000000E+00    7.920162281693801E-09
!   1.3500000E+02    8.689312433613942E-05
!                   -1.638765583241871E-10    0.000000000000000E+00   -5.445976435581539E-07
!                   -3.260546215008756E-09    0.000000000000000E+00    2.829427591902128E-12
!                   -3.351998343546536E-09    0.000000000000000E+00    6.612546761117837E-09
!   1.4500000E+02    8.175942046629851E-05
!                   -1.384428279847095E-10    0.000000000000000E+00   -4.838801342575158E-07
!                   -2.750557418695443E-09    0.000000000000000E+00    2.280796792454220E-12
!                   -2.818773289387068E-09    0.000000000000000E+00    5.569327072645579E-09
!   1.5500000E+02    7.718433330381285E-05
!                   -1.178576263805883E-10    0.000000000000000E+00   -4.325366577204661E-07
!                   -2.338664157789807E-09    0.000000000000000E+00    1.854285025140093E-12
!                   -2.390291633013819E-09    0.000000000000000E+00    4.728959785590588E-09
!   1.6500000E+02    7.308344678779011E-05
!                   -1.010541140283743E-10    0.000000000000000E+00   -3.887784140609058E-07
!                   -2.003069219405086E-09    0.000000000000000E+00    1.519812238302657E-12
!                   -2.042663757693794E-09    0.000000000000000E+00    4.045729979582435E-09
!   1.7500000E+02    6.938815426986298E-05
!                   -8.722588597628442E-11    0.000000000000000E+00   -3.512132394062391E-07
!                   -1.727335616477760E-09    0.000000000000000E+00    1.255906451289698E-12
!                   -1.758079035653547E-09    0.000000000000000E+00    3.485411979109969E-09
!   1.8500000E+02    6.604220927845390E-05
!                   -7.575765865429726E-11    0.000000000000000E+00   -3.187473619604029E-07
!                   -1.498995215177525E-09    0.000000000000000E+00    1.045428043545681E-12
!                   -1.523136694338215E-09    0.000000000000000E+00    3.022127700543810E-09
!   1.9500000E+02    6.299911914505276E-05
!                   -6.617793643013534E-11    0.000000000000000E+00   -2.905135869657135E-07
!                   -1.308481303347506E-09    0.000000000000000E+00    8.765125655926792E-13
!                   -1.327638445985104E-09    0.000000000000000E+00    2.636124214433246E-09
!   2.0500000E+02    6.022015719417070E-05
!                   -5.812005495505718E-11    0.000000000000000E+00   -2.658183543989876E-07
!                   -1.148424994398436E-09    0.000000000000000E+00    7.396755563497527E-13
!                   -1.163772150947790E-09    0.000000000000000E+00    2.312197498129573E-09
!   2.1500000E+02    5.767283276742131E-05
!                   -5.129894108907559E-11    0.000000000000000E+00   -2.441023741318757E-07
!                   -1.013065908463235E-09    0.000000000000000E+00    6.283734005544198E-13
!                   -1.025472579591845E-09    0.000000000000000E+00    2.038536126255646E-09
!   2.2500000E+02    5.532970311618965E-05
!                   -4.549044064644129E-11    0.000000000000000E+00   -2.249110992991893E-07
!                   -8.978890411845543E-10    0.000000000000000E+00    5.364593306069736E-13
!                   -9.080018056729546E-10    0.000000000000000E+00    1.805888654870618E-09
!   2.3500000E+02    5.316744283951857E-05
!                   -4.051510235621125E-11    0.000000000000000E+00   -2.078723716149099E-07
!                   -7.993187908752558E-10    0.000000000000000E+00    4.610017534929314E-13
!                   -8.076273023942558E-10    0.000000000000000E+00    1.606946455241786E-09
!   2.4500000E+02    5.116610904202790E-05
!                   -3.623089786692503E-11    0.000000000000000E+00   -1.926793662914099E-07
!                   -7.145049240937218E-10    0.000000000000000E+00    3.979006828440654E-13
!                   -7.213785011799856E-10    0.000000000000000E+00    1.435880398535937E-09
!   2.5500000E+02    4.930855649651253E-05
!                   -3.252367186377769E-11    0.000000000000000E+00   -1.790774573945117E-07
!                   -6.411533743202616E-10    0.000000000000000E+00    3.452533832472083E-13
!                   -6.468769976923290E-10    0.000000000000000E+00    1.288030321855556E-09
!   2.6500000E+02    4.757996872078616E-05
!                   -2.929998531591678E-11    0.000000000000000E+00   -1.668540276783255E-07
!                   -5.774118135315439E-10    0.000000000000000E+00    3.006504582226094E-13
!                   -5.822096625358964E-10    0.000000000000000E+00    1.159621502090200E-09
!   2.7500000E+02    4.596747933454624E-05
!                   -2.648481734175485E-11    0.000000000000000E+00   -1.558305006477632E-07
!                   -5.217721829911791E-10    0.000000000000000E+00    2.631475967484596E-13
!                   -5.258162497034863E-10    0.000000000000000E+00    1.047588811765136E-09
!   2.8500000E+02    4.445986426369507E-05
!                   -2.401586769430745E-11    0.000000000000000E+00   -1.458560633796081E-07
!                   -4.729990065137357E-10    0.000000000000000E+00    2.313102070370440E-13
!                   -4.764254683839667E-10    0.000000000000000E+00    9.494264761434548E-10
!   2.9500000E+02    4.304728994459756E-05
!                   -2.184207958790653E-11    0.000000000000000E+00   -1.368026879008472E-07
!                   -4.300722232546555E-10    0.000000000000000E+00    2.039687899127268E-13
!                   -4.329933246963261E-10    0.000000000000000E+00    8.630634503472928E-10
!  The total time is   3.56755400000000      seconds.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!