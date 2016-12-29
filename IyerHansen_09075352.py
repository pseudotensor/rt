import scipy,numpy
from scipy import *
from numpy import *
from sympy.mpmath import ellippi as Pi

# KNOWN PROBLEMS:
#    - BREAKS DOWN FOR a!=0 Kerr :-S
#      -> Look at Omega_minus:
#         The expression given does not go to zero for a=0 as expected.
#         Setting Omega_minus=0 gives correct answer for a=0 and 
#         also more reasonable answers for a>~0

def deflection_angle_Iyer2009(a,b,VERBOSE=False):

    '''Given BH spin (a) and impact parameter (b), compute deflection angle of null geodesics.
       REFERENCE: https://arxiv.org/pdf/0907.5352v1.pdf
       SEE ALSO: https://arxiv.org/pdf/0908.0085v1.pdf

       WIP: TRACKING OF SIGN FOR PRO-/RETROGRAD NOT ELEGANT YET...

    '''

    # Iyer+2009 say they track retro/prograde orbit via the sign of b
    a=sign(b)*a # necessary for retrograd (need to understand that better)
    
    SANITY=True

    s = sign(b*a) # direct orbits: s=+1 retrograd orbits: s=-1
    if a==0:
        s=1 # sign(b)
    b_s = s*abs(b) # b_s: positive magnitude of impact parameter
    b=b_s # this gives k^2=k^2_v2
    b_sc = -a + s*6*cos(arccos(-s*a)/3.)

    if (b<b_sc and b>0) or (b>b_sc and b<0):
        return nan # captured orbit

    if VERBOSE:
        print "Computing deflection angle for (b,b',a)=(",b,1.-s*b_sc/b,a,")"

    r_0 = 2.*b/sqrt(3.)*sqrt(1.-a**2/b**2)*cos(1./3.*arccos(-3.*sqrt(3.)*(1.-a/b_s)**2 /(b*(1.-(a/b)**2)**1.5) )) # eq (20)  in https://arxiv.org/pdf/0907.5352v1.pdf

    P = r_0*(1.+a/b_s)/(1.-a/b_s) # eq (16) in https://arxiv.org/pdf/0907.5352v1.pdf 
    Q = sqrt((P-2.)*(P+6.)) # eq (19) in https://arxiv.org/pdf/0907.5352v1.pdf & after eq (6) in https://arxiv.org/pdf/0611086v2.pdf

    omega_s = a/b_s # eq (36)
    h = 1./r_0 # eq (36)
    h_sc = (1.+omega_s)/(1.-omega_s) # eq (37) (~> 1 for a=0)
    omega_0 = a*a # eq (36)

    r_0_over_Q_v2 = 1. / h_sc / sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (40)  
    r_0_over_Q = r_0/Q # NOT THE SAME AS r_0_over_Q = 1./h_sc # eq(40)

    # SANITY CHECK    
    if SANITY and not isclose(r_0_over_Q,r_0_over_Q_v2):
        print "r0/Q=",r_0_over_Q,"r0/Q_v2",r_0_over_Q_v2," should equal each other" # see e.g. page 7 after eq (20)
    r_0_over_Q=r_0_over_Q_v2
    if VERBOSE:
        print "Using r_0__over_Q_v2"

    k_squared = (sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))+6.*h/h_sc-1.)/2./sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (41)
    k_squared_v2 = (Q-P+6.)/2./Q # see between eq(34) and eq(35) 
    
    # SANITY CHECK
    if SANITY and not isclose(k_squared,k_squared_v2):
        print "k^2=",k_squared,"k^2_v2",k_squared_v2," should equal each other and be within [0,1]" # see page 11 after eq (34)

    k_squared=k_squared_v2
    k=sqrt(k_squared)

    psi = arcsin(sqrt( ( 1.- 2.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) /  ( 1.- 6.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) ) # eq (42)
    psi_v2 = arcsin( sqrt( (Q+2.-P)/(Q+6.-P) ) ) # p11 before eq (35)

    # SANITY CHECK
    if SANITY and not isclose(psi,psi_v2):
        print "psi=",psi,"psi_v2",psi_v2," should equal each other" # see page 11 before eq (35)
    # psi=psi_v2; print "USING psi_v2 which differs from psi"

    # Omega_minus expression is wrong ...
    Omega_plus  =  +((1.+sqrt(1.-omega_0))*(1.-omega_s)-omega_0/2.)/sqrt(1.-omega_0)/(1.+sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)
    Omega_minus =  -((1.-sqrt(1.-omega_0))*(1.-omega_s)+omega_0/2.)/sqrt(1.-omega_0)/(1.-sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)

    u1 = (P-2.-Q)/4./r_0 # eq (13)
    u2 = 1/r_0 # eq (14)
    u_plus = (1.+sqrt(1.-a**2))/a**2 # eq (30)
    u_minus = (1.-sqrt(1.-a**2))/a**2 # eq (30)
    C_plus = ( (1.+sqrt(1.-a**2))*2.*(1.-omega_s)-a**2 ) / (2.*a**2*sqrt(1.-a**2)) # eq 32
    C_minus = ( a**2-(1.-sqrt(1.-a**2))*2.*(1.-omega_s) ) / (2.*a**2*sqrt(1.-a**2)) # eq 33

    Omega_plus_v2 = C_plus/(u_plus-u1) # see eq (32), eq after (34) at the bottom of page 11
    Omega_minus_v2 = C_minus/(u_minus-u1) # see eq (33), eq after (34) at the bottom of page 11

    if SANITY and not isclose(Omega_plus,Omega_plus_v2):
    # if Omega_plus!=Omega_plus_v2:
        print "Omega_plus=",Omega_plus,"Omega_plus_v2",Omega_plus_v2," should equal each other"
    if SANITY and not isclose(Omega_minus,Omega_minus_v2):
    # if Omega_minus!=Omega_minus_v2:
        KNOWN_PROBLEM=True
        if not KNOWN_PROBLEM:
            print "Omega_minus=",Omega_minus,"Omega_minus_v2",Omega_minus_v2," should equal each other"
            print "USING Omega_minus_v2"
        Omega_minus=Omega_minus_v2 # ; print "USING Omega_minus_v2"
    # C'mon THE DENOMINATOR IS ZERO FOR SCHWARZSCHILD!!!
    # print 1.-sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) # eq (43)

    try:
        n_plus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.+sqrt(1.-omega_0))) # eq (44)
        n_minus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.-sqrt(1.-omega_0))) # eq (44)
    except:
        n_plus=0.
        n_minus=0.
        pass

    n_plus_v2 = (u2-u1)/(u_plus-u1) # bottom of page 11 after eq (34)
    n_minus_v2 = (u2-u1)/(u_minus-u1) # bottom of page 11 after eq (34)

    # some more SANITY CHECKS
    if SANITY and not isclose(n_plus,n_plus_v2):
        print "n_plus=",n_plus,"n_plus_v2",n_plus_v2," should equal each other"
    if SANITY and not isclose(n_minus,n_minus_v2):
        print "n_minus=",n_minus,"n_minus_v2",n_minus_v2," should equal each other"

    if a==0.:
        if SANITY and Omega_minus!=0:
            print "Attention Omega_minus!=0 but a=0. Setting it to 0 now..."
            Omega_minus = 0. # maintext before eq (35)
        if SANITY and Omega_plus!=1:
            print "Attention Omega_plus!=1 but a=0. Setting it to 1 now..."
            Omega_plus = 1. # maintext before eq (35)
        if SANITY and n_minus!=0:
            print "Attention n_minus!=0 but a=0. Setting it to 0 now..."
            n_minus = 0. # maintext before eq (35) # no mention but given that n_plus=0 n_minus should be as well
        if SANITY and n_plus!=0:
            print "Attention n_plus!=0 but a=0. Setting it to 0 now..."
            n_plus = 0. # maintext before eq (35)

    from scipy.special import ellipk,ellipkinc

    deflection_angle_09075352v1_schwarzschild = -pi+4.*sqrt(r_0_over_Q) * ( ellipk(k_squared) - ellipkinc(psi,k_squared) ) # eq(35) in arXiv:0907.5352v1 

    deflection_angle_09075352v1 = -pi+4./(1.-omega_s)*sqrt(r_0_over_Q)*(Omega_plus*(Pi(n_plus,k_squared)-Pi(n_plus,psi,k_squared))+Omega_minus*(Pi(n_minus,k_squared)-Pi(n_minus,psi,k_squared)))   # eq (34) in https://arxiv.org/pdf/0907.5352v1.pdf

    if deflection_angle_09075352v1.imag!=0.:
        print "WARNING EXACT DEFLECTION ANGLE COMPLEX"
        deflection_angle_09075352v1 = float(deflection_angle_09075352v1.real)
    else:
        deflection_angle_09075352v1 = float(deflection_angle_09075352v1)

    deflection_angle = deflection_angle_09075352v1
    
    return deflection_angle
