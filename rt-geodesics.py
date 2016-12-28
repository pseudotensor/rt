## DIAGNOSTIC ROUTINE TO ANALYZE GEODESICS COMPUTED WITH ASTRORAY

# USAGE:
# python -i rt-geodesics.py [USER-FILENAME] 0 0
# [USER-FILENAME]: geodesics93.75th157fn5550geodesic*.dat
# Last two numbers select range (min,max) in viewing angles

# matplotlib.use('GTKAgg')
# matplotlib.use('Agg')
import string,matplotlib,glob,sys

try:
    __IPYTHON__ ## ARE WE RUNNING FROM IPYTHON?
except:
    matplotlib.use("Agg") # produces png without X
    # ioff() UNKNOWN WITHOUT IPYTHON

# from matplotlib.pyplot import *
# from numpy import * # loadtxt,savetxt,amin
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rc('font', size=15)

angle_min=int(sys.argv[-2])
angle_max=int(sys.argv[-1])
sys.argv = sys.argv[:-2] # discard two last int arg for viewing angle


## USER SPECS ##
HOME="/home/rgold/"
anim=False # figure(2)
figure1=False
figure4=False
figure3=False
figure5=True # deflection angle vs impact parameter


def b_sc(a,s): 
    '''Given black hole spin a and sign for pro/retrograde orbit s return b_sc as given in eq (28) from https://arxiv.org/pdf/0907.5352v1.pdf'''
    return -a+s*6.*cos(1./3.*arccos(-s*a))


def xyz2deflection(x,y,z,every=5): # CAREFUL! RESULTS DEPEND SOMEWHAT ON VALUE every ...
    '''IN: x,y,z coordinates along a geodesic
       OUT: deflection angle'''
    tangent_in  = array((diff(x[::every])[ 0],diff(y[::every])[ 0],diff(z[::every])[ 0]))
    tangent_out = array((diff(x[::every])[-1],diff(y[::every])[-1],diff(z[::every])[-1]))
    # A.B=|A||B|cos(angle)
    cos_angle = dot(tangent_in,tangent_out)/norm(tangent_in)/norm(tangent_out)
    return arccos(cos_angle)


def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    #u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j] #TAKES VERY LONG
    x=outer(np.cos(u),np.sin(v))
    y=outer(np.sin(u),np.sin(v))
    z=outer(ones(size(u)),np.cos(v))
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)


# FORMAT: ppy[geodesic_label].cooxx[p][geo_idx]
#          r=y[4]
# 0.000000 500.209461 0.895617 0.034918 500.000000 -499.790539 0.024889 0.032990 1 0
# 0.234413 499.975146 0.895612 0.034934 500.000941 -499.790341 0.024912 0.033020 1 1
#     ^^^ = geodesic_label geo_idx

# y[0] = tin;    
# y[1]= phiin; 
# y[2]=rprin; 
# y[3]=muprin; 
# y[4]=rin; # y[4] = r
# y[5]=muin;# y[5] = costh 
# y[6]=tprin;
# y[7]=phprin;
# y[8]=tain;
# y[9]=rain;
# y[10]=muain;
# y[11]=pain;

if figure1:
    fig1 = figure(1)
    ax1 = fig1.add_subplot(111)
if figure3:
    fig3 = figure(3)
    ax3 = fig3.add_subplot(111)
if figure4:
    fig4 = figure(4)
    ax4 = fig4.add_subplot(111)

fig3d = figure(2)
ax3d = fig3d.gca(projection='3d')

veto = 0 # counter for vetos on bad geodesics
bmax = 0 # max impact parameter of all geodesics
alpha_vs_b = []


## GET DATA ##

b_vs_deflection_angle_vrt2 = []

# GET GRTRANS DATA (SEE LATER: grtrans_geodesics_*.txt )

# d_grtrans = loadtxt("grtrans_alpha_vs_b.dat",comments='#')
# d_grtrans = d_grtrans[d_grtrans[:,4]<b_sc(a,-1)] # filter out captured orbits
# sign_of_b_grtrans = d_grtrans[:,2] d_grtrans[:,3] 
# d_eq_grtrans = d_grtrans[abs(d_grtrans[:,5])==amin(abs(d_grtrans[:,5])),:] # FOCUS ON SMALLEST Z COORDINATE: EQUATOR?
# d_eq_grtrans = d_grtrans[-20:,:]

# GET VRT2 DATA
d_vrt2=[]
# GET ODYSSEY DATA (Developer Hung-Yi Pu)
d_odyssey=[]
deflection_angle_odyssey=[]

# d_odyssey = loadtxt(HOME+"codes/odyssey/geodesic_output_alpha_x_y_z.txt")
for ODYSSEY_FILE in glob.glob(HOME+"codes/odyssey/files/data_*.txt"):
    d_odyssey.append(loadtxt(ODYSSEY_FILE))

    # d_odyssey = array(d_odyssey)
    if abs(d_odyssey[-1][0,0])>abs(b_sc(0.9375,-sign(d_odyssey[-1][0,0]))):
        # print "ODYSSEY CAPTURE ORBIT CHECK: ",d_odyssey[-1][0,0],b_sc(0.9375,sign(d_odyssey[-1][0,0]))
        deflection_angle_odyssey.append((d_odyssey[-1][0,0],xyz2deflection(d_odyssey[-1][:,1],d_odyssey[-1][:,2],d_odyssey[-1][:,3])))

deflection_angle_odyssey=array(deflection_angle_odyssey)

# WIP: COMPUTE DEFLECTION ANGLE FROM x,y,z and compare to grtrans_alpha_vs_b.dat 9th col
# print "d_grtrans[-1,8]=",d_grtrans[-1,8]," deflection_angle=",xyz2deflection(d_grtrans[-1,2],d_grtrans[-1,3],d_grtrans[-1,4])


# GET vrt data (Avery's code) lambda t x y z
for FILE_vrt2 in glob.glob(HOME+"codes/grtrans/avery_ray_tracing_test_data/eta_*_n.d"):
    d_vrt2.append(loadtxt(FILE_vrt2,skiprows=1,usecols=(2,3,4))) # get x,y,z coordinates of ray trajectory
    a=0.9375 ## HARDWIRED
    x_vrt2,y_vrt2,z_vrt2 = d_vrt2[-1][:,0],d_vrt2[-1][:,1],d_vrt2[-1][:,2]
    # if sqrt(x_vrt2[-1]**2+y_vrt2[-1]**2)<b_sc(a,1):
    if amin(sqrt(x_vrt2**2+y_vrt2**2))<b_sc(a,1):
        # y_vrt2[0]<b_sc(a,1):
        continue # pass
    b_vrt2 = y_vrt2[0] # impact parameter

    deflection_angle_vrt2 = xyz2deflection(x_vrt2,y_vrt2,z_vrt2)
    # print "Deflection angle (vrt2):",deflection_angle_vrt2,"rad =",deflection_angle_vrt2/2./pi*360.,"deg b=",b_vrt2
    
    # b_vs_cos_angle_vrt2 = hstack((b_vs_cos_angle_vrt2,(cos_angle_vrt2,b_vrt2)))
    b_vs_deflection_angle_vrt2.append((b_vrt2,deflection_angle_vrt2))
    
b_vs_deflection_angle_vrt2 = array(b_vs_deflection_angle_vrt2)
alpha_vs_b_darwin = []
alpha_vs_b_09075352 = []

every_nth_geodesic=1
every_nth_point_on_each_geodesic=1 # FIXME: UNEXPECTED BEHAVIOR FOR VALUES != 1
for FILE in sys.argv[1:][::every_nth_geodesic]:
    filename=FILE # e.g. ~> 93.75th43fn2150geodesic.dat

    d=loadtxt(filename)[::every_nth_point_on_each_geodesic,:]

    a=float(filename.split("/")[-1].split("geodesics")[1].split("th")[0])/100. # spin from filename
    # a=0.001; print "[HARDWIRED]: a=",a # deflection angles offset by ~3.36rad independent of b
    
    rh=1.+sqrt(1.-a*a);

    lambda_affine=d[:,8]

    # ASTRORAY [geoint.cpp]: y= t, phi, rp ,mup, r, mu ,tpr,phpr
    # outputting ppy[currth].cooxx[p][geo_idx] to file
    
    # CORRECT: see evalpointzero: cooxx <-> ly
    r = d[:,1]
    costheta = d[:,2]
    phi = d[:,3]

    # r[find(r<=rh)[0]:]=None # WHAT HAPPENS INSIDE A BH STAYS INSIDE A BH
    
    # chuck off geodesics that entered BH
    HIDE_GEODESICS_INSIDE_HORIZON=True
    if HIDE_GEODESICS_INSIDE_HORIZON==True:
        entered_horizon_mask = find( r <= rh )
    else:
        entered_horizon_mask = find( r == r ) # Always True
    # how to merge the masks?

    r        = delete(r,entered_horizon_mask)    
    costheta = delete(costheta,entered_horizon_mask)    
    phi      = delete(phi,entered_horizon_mask)    

    lambda_affine = delete(lambda_affine,entered_horizon_mask) # r,th,ph arrays have changed
    runaway_mask = find(lambda_affine > 1.5)
    r        = delete(r,runaway_mask)    
    costheta = delete(costheta,runaway_mask)    
    phi      = delete(phi,runaway_mask)    
    lambda_affine = delete(lambda_affine,runaway_mask)

    # last elements seem off...
    # r=r[:-1]
    # costheta=costheta[:-1]
    # phi=phi[:-1]

    # print "sizes (b4):",size(r),size(HORIZON_MASK)
    # r2=r[HORIZON_MASK]
    # print "sizes (after):",size(r),size(HORIZON_MASK)

    # costheta=costheta[HORIZON_MASK]
    # phi=phi[HORIZON_MASK]

    # r[find(lambda_affine>=1.5)[0]:]=None # NOT INTERESTED IN FOLLOWING REMOTE GEODESICS

    sintheta = sqrt(1.-costheta*costheta)
    
    x = r*sintheta*cos(phi)
    y = r*sintheta*sin(phi)
    z = r*costheta

    # A.B=|A||B|cos(angle)
    deflection_angle_astroray = xyz2deflection(x,y,z)
    # print "Deflection angle ",deflection_angle_astroray,"rad =",deflection_angle_astroray/2./pi*360.,"deg"

    # http://arxiv.org/pdf/gr-qc/9907034v1.pdf eqs.(20),(24)
    # http://arxiv.org/pdf/1405.2919.pdf see Figs. 2.2,3.5
    # https://arxiv.org/pdf/0907.5352v1.pdf eq(34) exact bending angle in equatorial Kerr

    ##############################################################################################
    # !!! ===>   https://arxiv.org/pdf/gr-qc/0611086v2.pdf eqs.(8),(14->23),(17->25)    <=== !!! #
    ##############################################################################################
    r_per=amin(r) # For Schwarzschild b=r_per/sqrt(1.-r_s/r_per) http://arxiv.org/abs/gr-qc/9907034
    b=d[-1,12]

    # b_crit = 3.*sqrt(3.) # "photon sphere":3 https://arxiv.org/pdf/gr-qc/0611086v2.pdf eq 

    s = sign(b*a) # direct orbits: s=+1 retrograd orbits: s=-1
    b_s = s*abs(b) # b_s: positive magnitude of impact parameter
    b=abs(b); # print "[WARNING]: ENFORCE b>=0" # should really use b_s ...
    
    # b_sc=-a+s*6.*cos(1./3.*arccos(-s*a))
    b_crit = b_sc(a,s)
    # b_prime = 1.-s*b_sc/b 
    # b_prime = 1.-s*b_sc/b 
    b_prime = 1.-b_crit/abs(b) # https://arxiv.org/pdf/gr-qc/0611086v2.pdf eq before eq (20) III.D page 12 ? # abs(b)?

    # elliptic integrals: http://arxiv.org/pdf/math/9409227v1.pdf
    deflection_angle_0611086v2  = 0. # Schwarzschild need elliptic integral of the 1st kind

    try:
        from sympy.mpmath import ellippi as Pi
        EXACT_BENDING_ANGLE=True
    except ImportError:
        EXACT_BENDING_ANGLE=False
        continue
    # EXACT_BENDING_ANGLE=False # ...WIP... NOT WORKING AS EXPECTED YET

    sys.path.append("/home/rgold/codes/astroray/")
    from IyerHansen_09075352 import deflection_angle_Iyer2009

    # https://arxiv.org/pdf/0907.5352v1.pdf

    # r_0 = r_per # see eq(20) in https://arxiv.org/pdf/0907.5352v1.pdf
    r_0 = 2.*abs(b)/sqrt(3.)*sqrt(1.-a**2/abs(b)**2)*cos(1./3.*arccos(-3.*sqrt(3.)*(1.-a/b_s)**2 /(b*(1.-(a/abs(b))**2)**1.5) )) # eq (20)  in https://arxiv.org/pdf/0907.5352v1.pdf

    # What about b_s=0 ?  :-s
    P = r_0*(1.+a/b_s)/(1.-a/b_s) # eq (16) in https://arxiv.org/pdf/0907.5352v1.pdf 
    Q = sqrt((P-2.)*(P+6.)) # eq (19) in https://arxiv.org/pdf/0907.5352v1.pdf

    omega_s = a/b_s # eq (36)
    h = 1./r_0 # eq (36)
    h_sc = (1.+omega_s)/(1.-omega_s) # eq (37) (~> 1 for a=0)
    omega_0 = a*a # eq (36)

    # r_0_over_Q = 1. / h_sc / sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (40)  
    r_0_over_Q = r_0/Q 

    k_squared = (sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))+6.*h/h_sc-1.)/2./sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (41)
    k_squared_v2 = (Q-P+6.)/2./Q # see between eq(34) and eq(35) DOES NOT AGREE with eq(41) k_squared above!
    # print "k_squared,k_squared_v2:",k_squared,k_squared_v2
    k_squared=k_squared_v2
    k=sqrt(k_squared)
    psi = sqrt( ( 1.- 2.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) /  ( 1.- 6.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) )) # eq (42)
    Omega_plus  =  +((1.+sqrt(1.-omega_0))*(1.-omega_s)-omega_0/2.)/sqrt(1.-omega_0)/(1.+sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)
    Omega_minus  =  -((1.-sqrt(1.-omega_0))*(1.-omega_s)+omega_0/2.)/sqrt(1.-omega_0)/(1.-sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)
    n_plus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.+sqrt(1.-omega_0))) # eq (43)
    n_minus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.-sqrt(1.-omega_0))) # eq (43)

    from scipy.special import ellipk,ellipkinc
    # FIXME: SOMETHING IS NOT RIGHT HERE, THIS SHOULD AT LEAST AGREE WITH darwin_tab[:,-1] ...
    deflection_angle_09075352v1_schwarzschild = -pi+4.*sqrt(r_0_over_Q) * ( ellipk(k_squared) - ellipkinc(psi,k_squared) ) # eq(35) in arXiv:0907.5352v1 
    alpha_vs_b_09075352.append([s*abs(b),deflection_angle_09075352v1_schwarzschild])

    if EXACT_BENDING_ANGLE==True:
        # try:
            from KerrDeflection import SchwarzDeflection,FindRoots,EquatorialDeflection # ,SchwTiltCoords
            b_exact_grid = arange(5,21)
            # b_exact_grid = arange(-20,-10) # DOES NOT WORK, HAVE TO PLAY WITH INCLINATION INSTEAD
            # b_exact_grid = array([10])
            bignumber=1e20
            deflection_angle_schwarzschild_exact = SchwarzDeflection(bignumber,b_exact_grid)-pi
            # deflection_angle_schwarzschild_exact = SchwTiltCoords(deflection_angle_schwarzschild_exact+pi,pi,b_exact_grid,0)[1]

            deflection_angle_kerr_eq_exact = []
            for b_tmp in b_exact_grid: 
                roots=FindRoots(a,pi,bignumber,array([b_tmp]),0) # spin,i,E,b,?
                deflection_angle_kerr_eq_exact.append(EquatorialDeflection(a, bignumber, b_tmp, roots.reshape(4,-1))[0]-pi)
            # roots=FindRoots(0,pi,bignumber,b_exact_grid,0) # spin,i,E,b,?
            # deflection_angle_kerr_eq_exact = EquatorialDeflection(0, bignumber, b_exact_grid, roots.reshape(4,-1))

            # UNDERSTANDING CONVENTIONS AND FCT ARGS: THE FOLLOWING PRODUCES CONSISTENT RESULTS
            # KerrDeflection(0,pi,bignumber,arange(21),0.)[1],SchwarzDeflection(bignumber,arange(21))-pi

        # except:
        #     pass

    if EXACT_BENDING_ANGLE==True and k==k:
        # k^2
        deflection_angle_09075352v1 = -pi+4./(1.-omega_s)*sqrt(r_0_over_Q)*(Omega_plus*(Pi(n_plus,k_squared)-Pi(n_plus,psi,k_squared))+Omega_minus*(Pi(n_minus,k_squared)-Pi(n_minus,psi,k_squared)))   # eq (34) in https://arxiv.org/pdf/0907.5352v1.pdf
        # deflection_angle_09075352v1 = -pi+4./(1.-omega_s)*sqrt(r_0_over_Q)*(Omega_plus*(Pi(n_plus,k)-Pi(n_plus,psi,k))+Omega_minus*(Pi(n_minus,k)-Pi(n_minus,psi,k)))   # eq (34) in https://arxiv.org/pdf/0907.5352v1.pdf
        if deflection_angle_09075352v1.imag!=0.:
            print "WARNING EXACT DEFLECTION ANGLE COMPLEX"
            deflection_angle_09075352v1 = float(deflection_angle_09075352v1.real)
        else:
            deflection_angle_09075352v1 = float(deflection_angle_09075352v1)

    deflection_angle_weakfield   = 4./b * ( 1. + 15./16.*pi/b )
    deflection_angle_weakfield_Sereno = 4./b + (15.*pi/4.-4.*a)/b**2 + (4*a**2-10.*pi*a+128./3.)/b**3 + (15./64.*pi*(76*a**2+231)-4.*a*(a**2+48.))/b**4 + (4.*(a**2+128)*a**2-9./2.*pi*(6.*a**2+77.)*a+3584./5.)/b**5 # http://arxiv.org/pdf/1405.2919.pdf eq (3.40) (equatorial), see: Sereno,de Luca for general orbits

    deflection_angle_strongfield = log( 3.482/(b-3.*sqrt(3.))) # http://arxiv.org/pdf/gr-qc/9907034v1.pdf SCHWARZSCHILD
    deflection_angle_0611086v2 = -pi + log( 216.*(7.-4.*sqrt(3.))/b_prime) + (-17.+4.*sqrt(3.)+5.*log(216.*(7.-4.*sqrt(3.))/b_prime))*b_prime/18. + (-879.+236.*sqrt(3.)+205.*log(216.*(7.-4.*sqrt(3.))/b_prime))*b_prime**2/1296. + (-321590.+90588*sqrt(3.)+68145.*log(216.*(7.-4.*sqrt(3.))/b_prime))*b_prime**3/629856. # SCHWARZSCHILD
    #    try:
    if True:
        # x_0 = b # P in Darwin is perihelion distance...
        x_0 = r_per # P in Darwin is perihelion distance...
        darwin_lambda = (3.- x_0 - sqrt(x_0**2+2*x_0-3.))/(3.- x_0 + sqrt(x_0**2+2*x_0-3.))
        darwin_phi_0 = sqrt( (-3.+x_0-sqrt(x_0**2+2*x_0-3.) ) / (2.*(2.*x_0-3.)) )
        darwin_G = sqrt( ( 8.*x_0*(-3.+x_0+sqrt(x_0**2+2.*x_0-3.)) ) / (2.*x_0-3.) )
        # 2DO: F<->python fct incomplete elliptic integral of the first kind
        from scipy.special import ellipkinc
        F_darwin=ellipkinc(darwin_phi_0,darwin_lambda) # darwin_lambda <-> k wiki <-> m in scipy
        deflection_angle_darwin_exact = -pi + F_darwin*darwin_G # Darwin 1959
        # print "darwin_lambda,darwin_phi_0,darwin_G,x_0,F_darwin,deflection_angle_darwin_exact:",darwin_lambda,darwin_phi_0,darwin_G,x_0,F_darwin,deflection_angle_darwin_exact
        deflection_angle_darwin = -pi + 2.*log( 36.*(2.-sqrt(3.))/(r_per-3.)) # Darwin 1959
        alpha_vs_b_darwin.append([s*abs(b),deflection_angle_darwin_exact])

#    except:
#        pass
    # print "pericenter=",r_per,"impact parameter b=",b
    # print "...assuming Schwarzschild for weakfield..."
    # print "deflection_angle_weakfield=",deflection_angle_weakfield,"rad =",deflection_angle_weakfield/2./pi*360.,"deg"
    # print "deflection_angle_weakfield_Sereno=",deflection_angle_weakfield_Sereno,"rad =",deflection_angle_weakfield_Sereno/2./pi*360.,"deg"
    # print "deflection_angle_strongfield=",deflection_angle_strongfield,"rad =",deflection_angle_strongfield/2./pi*360.,"deg"
    # print "deflection_angle_darwin=",deflection_angle_darwin,"rad =",deflection_angle_darwin/2./pi*360.,"deg"

    if diff(r)[-1]>0: # HYPERBOLIC ORBIT
        alpha_vs_b += [[s*abs(b),deflection_angle_astroray,deflection_angle_weakfield,deflection_angle_weakfield_Sereno,deflection_angle_09075352v1]] #deflection_angle_0611086v2]] # deflection_angle_darwin]]

    veto += size(find(r==None))
    if False: # (amax(r)>500 and size(lambda_affine[r==amax(r)])>1):
        veto += 1
    else:
        ax3d.plot(x,y,z,"r-",alpha=0.2,linewidth=1)

    ax3d.plot(x_vrt2,y_vrt2,z_vrt2,"b-",alpha=0.9,linewidth=1)


    if figure1:
        ax1.plot(lambda_affine,r)
        # ax1.plot(d[:,-1],d[:,4],'r+',label="5th col, geodesic "+str(int(d[0,-2])))

    if figure3==True:
        figure(3)
        ax3.plot(lambda_affine)

    if figure4==True:
        figure(4)
        ax4.plot(x,y)


    completeness=sys.argv[1:].index(FILE)/float(size(sys.argv[1:]))
    if (completeness*1000)%100==0:
        print "completeness =",int(completeness*100),"%"

print "THERE HAVE BEEN",veto,"VETOS."



fig6=figure(6)
ax_grtrans_2d=fig6.add_subplot(111)
grtrans_deflection_angle = []
grtrans_b=[]

GEODESICS=glob.glob("grtrans_geodesic_*_alpha*_beta*.txt")
list.sort(GEODESICS,key=lambda x : x.split('beta')[1].split('.txt')[0]) # sort according to BETAS
# GEODESICS=GEODESICS[-100:][::4] # WORKS GIVES AGREEMENT PLOT: GRTRANS,VRT2,ASTRORAY,ODYSSEY
GEODESICS=GEODESICS[-100:][::1]

for FILE in GEODESICS:
# for FILE in glob.glob("grtrans_geodesic_*_alpha1.3_beta*.txt"):
# for FILE in glob.glob("grtrans_geodesic_*_alpha*_beta1.3.txt"): # nearly equatorial
# for FILE in glob.glob("grtrans_geodesic_*_alpha1.3_beta1.3.txt"):

    alpha=float(FILE.split('alpha')[1].split('_')[0])
    beta=float(FILE.split('beta')[1].split('.txt')[0])
    x_grtrans,y_grtrans,z_grtrans = loadtxt(FILE)
    r_eq_grtrans = sqrt(alpha**2+beta**2)
    # print FILE,x_grtrans[0],y_grtrans[0],z_grtrans[0]

    # if amin(r_eq_grtrans)>b_sc(a,1): # filter out captured orbit

    if amin(r_eq_grtrans)>abs(b_sc(a,-sign(alpha))): # filter out captured orbit
    # if True:
        grtrans_b.append(alpha)
        grtrans_deflection_angle.append(xyz2deflection(x_grtrans,y_grtrans,z_grtrans))

    ax3d.plot(x_grtrans,y_grtrans,z_grtrans,"g-",alpha=0.5,linewidth=1)
    ax_grtrans_2d.plot(x_grtrans,y_grtrans,"g-",alpha=0.5,linewidth=1)
    ax_grtrans_2d.set_xlabel('X');ax_grtrans_2d.set_ylabel('Y')
    ax_grtrans_2d.axis('equal')

ax3d.plot(x,y,z,"r-",alpha=0.2,linewidth=1,label="ASTRORAY")
ax3d.plot(x_grtrans,y_grtrans,z_grtrans,"g-",alpha=0.5,linewidth=1,label="GRTRANS")
# ax3d.plot(d_grtrans[:,3],d_grtrans[:,4],d_grtrans[:,5],"g-",alpha=0.5,linewidth=1,label="GRTRANS")
ax3d.plot(x_vrt2,y_vrt2,z_vrt2,"b-",alpha=0.5,linewidth=1,label="VRT2")
ax3d.legend()
# LIM=100
# ax3d.set_xlim3d(-LIM,LIM)
# ax3d.set_ylim3d(-LIM,LIM)
# ax3d.set_zlim3d(-LIM,LIM)


(xs,ys,zs) = drawSphere(0.,0.,0.,rh)
ax3d.plot_surface(xs,ys,zs,color='k')
ax3d.set_aspect("equal")
ax3d.set_xlabel("X")
ax3d.set_ylabel("Y")
ax3d.set_zlabel("Z")


#################
## ANNOTATIONS ##

if figure1==True:
    figure(1)
    # ax1.set_xlabel("Point on geodesic")
    ax1.set_xlabel("Affine parameter "+r"$\lambda$"+" on geodesic")
    ax1.set_ylabel("r")
    # ax1.legend()
    # ax1.axis((0,1e4,0,None))

if figure3==True:
        ax3.set_ylabel(r"$\lambda$")
        ax3.set_title("Do I look monotonic and smooth?")

# d_grtrans = d_grtrans[d_grtrans[:,5]>b_sc(a,1)] # filter out captured orbits
# d_grtrans = d_grtrans[d_grtrans[:,2]>b_sc(a,1)] # filter out captured orbits

from scipy import interpolate

# UNSORTED
grtrans_deflection_angle_int_fct = interpolate.interp1d(grtrans_b,grtrans_deflection_angle,kind=3, bounds_error=False, fill_value=None) #returns interpolating function
vrt2_deflection_angle_int_fct = interpolate.interp1d(b_vs_deflection_angle_vrt2[:,0],b_vs_deflection_angle_vrt2[:,1],kind=3, bounds_error=False, fill_value=None) #returns interpolating function

# BEFORE SPLINE WE NEED TO SORT b,alpha BUT FOR interp1d unsorted seems ok
# grtrans_sort_idx = argsort(grtrans_b)
# grtrans_b_sorted = array(grtrans_b)[grtrans_sort_idx]
# grtrans_deflection_angle_sorted = array(grtrans_deflection_angle)[grtrans_sort_idx]

# grtrans_deflection_angle_int_sorted_fct = interpolate.interp1d(grtrans_b_sorted,grtrans_deflection_angle_sorted,kind=3, bounds_error=False, fill_value=None) #returns interpolating function


if figure5==True:

    # figure(5)
    figure(5,figsize=(10,12))
    subplot(311)

    grid(True,alpha=0.5) # want shadow to cover grid, but zorder has no effect...
    alpha_vs_b = array(alpha_vs_b)
    # IGNORE GRTRANS DATA UNTIL SETUP IS CORRECT (-> EQUATORIAL GEODESICS)
    # plot(d_grtrans[:,4],d_grtrans[:,5]/pi,'gd',label="GRTRANS")
    # plot(d_grtrans[:,6],d_grtrans[:,8]/pi,'gd',label="GRTRANS",alpha=0.2)
    plot(-array(grtrans_b),array(grtrans_deflection_angle)/pi,'gd',label="GRTRANS",alpha=0.5)

    # plot(d_eq_grtrans[:,5],d_eq_grtrans[:,8]/pi,'gd',label="GRTRANS")
    plot(alpha_vs_b[:,0],alpha_vs_b[:,1]/pi,'b+',label="ASTRORAY")
    plot(-b_vs_deflection_angle_vrt2[:,0],b_vs_deflection_angle_vrt2[:,1]/pi,'mx',label="VRT2",alpha=0.5) # minus sign due to convention difference in setup

    # plot(-d_odyssey[:,0],deflection_angle_odyssey/pi,'ro',label="ODYSSEY",alpha=0.5) # minus sign
    plot(-deflection_angle_odyssey[:,0],deflection_angle_odyssey[:,1]/pi,'ro',label="ODYSSEY",alpha=0.5) # minus sign due to convention difference in setup

    ## CHECK python code in appendix of A.1.1 and cpp code in A.1.2. of https://arxiv.org/pdf/1405.2919.pdf

    image_center_idx = find(alpha_vs_b[:,0]>0.)[0] # AVOID PLOTTING ARTEFACTS ACROSS SHADOW
    plot(alpha_vs_b[:image_center_idx,0],abs(alpha_vs_b[:image_center_idx,2])/pi,'.',color='gray',label="weak field limit")
    plot(alpha_vs_b[:image_center_idx,0],abs(alpha_vs_b[:image_center_idx,3])/pi,'--',color='gray',label="weak field w spin")

    # plot(alpha_vs_b[:,0],abs(alpha_vs_b[:,2])/pi,'.',color='gray',label="weak field limit")
    # plot(alpha_vs_b[:,0],abs(alpha_vs_b[:,3])/pi,'--',color='gray',label="weak field w spin")

    # DARWIN (SCHWARZSCHILD) EXACT
    alpha_vs_b_darwin = array(alpha_vs_b_darwin)
    darwin_tab = array([[3.2,5.23,273.],[3.4,5.30,205.],[3.6,5.40,162.],[3.8,5.53,143.],[4,5.66,125.],[5,6.46,79.],[6,7.35,58.],[7,8.28,46.],[8,9.24,38.],[9,10.22,32.],[10,11.20,28.],[11,12.17,25.],[12,13.15,23.]]) # TABLE 3 in http://www.jstor.org/stable/pdf/100508.pdf # P:Perihelion, l:impact parameter mu:deflection angle
    darwin_tab[:,-1] *= pi/180.

    plot(darwin_tab[:,1],darwin_tab[:,2]/pi,'kd',lw=2,label="DARWIN (exact, tab)")
    plot(-darwin_tab[:,1],darwin_tab[:,2]/pi,'kd',lw=2) # Schwarzschild: symmetric under l->-l
    # try:
    #     # plot(b_exact_grid, deflection_angle_schwarzschild_exact/pi,'k^-',lw=2,label="Schwarzschild (exact)")
    #     # plot(b_exact_grid, array(deflection_angle_kerr_eq_exact)/pi,'kv-',lw=2,label="Kerr (exact)")
    # except:
    #     pass

    # angle_Iyer2009 = array([deflection_angle_Iyer2009(0,b_tmp) for b_tmp in b_exact_grid[3:]])
    b_exact_grid = concatenate((arange(-20,-6.9),arange(6,20.1)))
    b_exact_grid_idx = find(b_exact_grid>0)
    angle_Iyer2009 = array([deflection_angle_Iyer2009(a,b_tmp) for b_tmp in b_exact_grid])
    plot(b_exact_grid[b_exact_grid_idx], angle_Iyer2009[b_exact_grid_idx]/pi,'k-',label="exact,Iyer+2009")
    plot(b_exact_grid[:b_exact_grid_idx[0]-1], angle_Iyer2009[:b_exact_grid_idx[0]-1]/pi,'k-')

    # plot(alpha_vs_b_darwin[:,0],alpha_vs_b_darwin[:,1]/pi,'kv',lw=2,label="DARWIN (exact)")
    alpha_vs_b_09075352=array(alpha_vs_b_09075352)
    # SEEMS WRONG
    # plot(alpha_vs_b_09075352[:,0],alpha_vs_b_09075352[:,1]/pi,'k-.',lw=2,label="09075352 (a=0)")
    plot(alpha_vs_b[image_center_idx:,0],abs(alpha_vs_b[image_center_idx:,2])/pi,'.',color='gray')
    plot(alpha_vs_b[image_center_idx:,0],abs(alpha_vs_b[image_center_idx:,3])/pi,'--',color='gray')

    # plot(-b_vs_deflection_angle_vrt2[:,0],grtrans_deflection_angle_int_fct(b_vs_deflection_angle_vrt2[:,0])/pi,'.',color='k',label="SPLINE TO GRTRANS") # JUST TO CHECK

    # plot(alpha_vs_b[:image_center_idx,0],alpha_vs_b[:image_center_idx,3]/pi,'k-',label="strong field limit")
    # plot(alpha_vs_b[image_center_idx:,0],alpha_vs_b[image_center_idx:,3]/pi,'k-')

    # sorted_b = alpha_vs_b[searchsorted(sort(alpha_vs_b[:,0]),alpha_vs_b[:,0]),0]
    # sorted_angles = alpha_vs_b[searchsorted(sort(alpha_vs_b[:,0]),alpha_vs_b[:,0]),2]
    # sorted_data = alpha_vs_b[searchsorted(sort(alpha_vs_b[:,0]),alpha_vs_b[:,0]),2]
    # plot(sort(alpha_vs_b[:,0]),sorted_data,'k--',label="weak field limit")
    # axvspan(-b_crit,b_crit,fc="k",alpha=0.5)
    axvspan(0-0.1,0+0.1,fc="k",alpha=0.5)
    axvspan(b_sc(a,-1),b_sc(a,1),fc="k",alpha=0.5)
    # text(-2,0.5,"shadow",color="k",rotation="vertical",fontsize=20)
    text(0.51,0.51,"shadow",color="k",rotation="vertical",fontsize=20,transform=gca().transAxes,verticalalignment='center')
    # # axis((-bmax,bmax,-0.05,amax(alpha_vs_b[:,1])/pi+0.05))
    axis((-20,20,-0.05,amax(alpha_vs_b[:,1])/pi+0.05))
    xticks(visible=False)
    # xlabel(u"impact parameter "+r"$/M$")
    ylabel(u"deflection angle "+r"$/\pi$")
    legend(labelspacing=0.1); setp(gca().get_legend().get_texts(), fontsize='12')


    ## FIT SPLINES TO DEFLECTION ANGLES SO WE CAN COMPUTE DIFFERENCES
    ## WORKS, BUT MIGHT WANT TO USE interp() instead so we only have to interpolate one data set

    # # import scipy.signal
    # from scipy import interpolate
    # # impact_parameter_grid = linspace(5,10,6) # 4debugging
    # # NOT VERY CONVENIENT, DON'T WANT TO INTERPOLATE ACROSS SHADOW
    # # impact_parameter_grid = concatenate((linspace(-10,b_sc(a,-1),10),linspace(b_sc(a,1),10,10)))
    # GET ORIENTATION CONVENTION RIGHT:
    # impact_parameter_grid_retrograd = linspace(-10,b_sc(a,-1),10)
    # impact_parameter_grid_prograd = linspace(b_sc(a,1),10,10)
    impact_parameter_grid_prograd = linspace(-10,-b_sc(a,1),10)
    impact_parameter_grid_retrograd = linspace(-b_sc(a,-1),10,10)

    # BEFORE SPLINE WE NEED TO SORT b,alpha
    grtrans_sort_idx = argsort(grtrans_b)
    grtrans_b_sorted = array(grtrans_b)[grtrans_sort_idx]
    grtrans_deflection_angle_sorted = array(grtrans_deflection_angle)[grtrans_sort_idx]
    vrt2_sort_idx = argsort(-b_vs_deflection_angle_vrt2[:,0])
    vrt2_b_sorted = array(-b_vs_deflection_angle_vrt2[:,0])[vrt2_sort_idx]
    vrt2_deflection_angle_sorted = array(-b_vs_deflection_angle_vrt2[:,1])[vrt2_sort_idx]

    grtrans_deflection_angle_retrograd_int = interpolate.splev(impact_parameter_grid_retrograd,interpolate.splrep(grtrans_b_sorted,grtrans_deflection_angle_sorted,k=5))
    grtrans_deflection_angle_prograd_int = interpolate.splev(impact_parameter_grid_prograd,interpolate.splrep(grtrans_b_sorted,grtrans_deflection_angle_sorted,k=5))

    vrt2_deflection_angle_retrograd_int = interpolate.splev(impact_parameter_grid_retrograd,interpolate.splrep(vrt2_b_sorted,vrt2_deflection_angle_sorted,k=5))
    vrt2_deflection_angle_prograd_int = interpolate.splev(impact_parameter_grid_prograd,interpolate.splrep(vrt2_b_sorted,vrt2_deflection_angle_sorted,k=5))

    # grtrans_deflection_angle_int_fct = interpolate.interp1d(grtrans_b,grtrans_deflection_angle,kind=3, bounds_error=False, fill_value=None) #returns interpolating function
    # grtrans_deflection_angle_prograd_int_fct = interpolate.interp1d(grtrans_b,grtrans_deflection_angle,kind=3) #returns interpolating function
    # grtrans_deflection_angle_retrograd_int_fct = interpolate.interp1d(grtrans_b,grtrans_deflection_angle,kind=3) #returns interpolating function

    # figure(50)

    # tight_layout(h_pad=0.1,w_pad=0)

    # savefig("deflection_angle_vs_impact_parameter.png")

    # exit()

    # figure(51)
    subplot(312)

    semilogy(alpha_vs_b[:,0],abs(vrt2_deflection_angle_int_fct(-alpha_vs_b[:,0])-alpha_vs_b[:,3])/alpha_vs_b[:,3],'cD',label="VRT2 vs WEAKFIELD")
    semilogy(-b_vs_deflection_angle_vrt2[:,0],abs(grtrans_deflection_angle_int_fct(b_vs_deflection_angle_vrt2[:,0])-b_vs_deflection_angle_vrt2[:,1])/b_vs_deflection_angle_vrt2[:,1],'g^',label="GRTRANS vs VRT2")
    semilogy(-deflection_angle_odyssey[:,0],abs(vrt2_deflection_angle_int_fct(deflection_angle_odyssey[:,0])-deflection_angle_odyssey[:,1])/deflection_angle_odyssey[:,1],'r.',label="ODYSSEY vs VRT2")
    semilogy(alpha_vs_b[:,0],abs(vrt2_deflection_angle_int_fct(-alpha_vs_b[:,0])-alpha_vs_b[:,1])/alpha_vs_b[:,1],'b+',label="ASTRORAY vs VRT2")
    # plot(impact_parameter_grid_retrograd,grtrans_deflection_angle_retrograd_int/pi,'k-',label="SPLINE")
    # plot(impact_parameter_grid_prograd,grtrans_deflection_angle_prograd_int/pi,'k-')
    axvspan(0-0.1,0+0.1,fc="k",alpha=0.5)
    axvspan(b_sc(a,-1),b_sc(a,1),fc="k",alpha=0.5)
    # text(0,0.1,"shadow",color="k",rotation="vertical",fontsize=20)
    text(0.51,0.51,"shadow",color="k",rotation="vertical",fontsize=20,transform=gca().transAxes,verticalalignment='center')
    grid(True)
    legend(labelspacing=0.1); setp(gca().get_legend().get_texts(), fontsize='12')
    axis((-20,20,None,None))
    # xlabel(u"impact parameter "+r"$/M$")
    xticks(visible=False)
    ylabel(r"$\Delta$"u" deflection angle (relative error)")


    subplot(313)
    semilogy(alpha_vs_b[:,0],abs(vrt2_deflection_angle_int_fct(-alpha_vs_b[:,0])-alpha_vs_b[:,3]),'cD',label="VRT2 vs Sereno+ (2014)") # weakfield with spin
    semilogy(-b_vs_deflection_angle_vrt2[:,0],abs(grtrans_deflection_angle_int_fct(b_vs_deflection_angle_vrt2[:,0])-b_vs_deflection_angle_vrt2[:,1]),'g^',label="GRTRANS vs VRT2")
    semilogy(-deflection_angle_odyssey[:,0],abs(vrt2_deflection_angle_int_fct(deflection_angle_odyssey[:,0])-deflection_angle_odyssey[:,1]),'r.',label="ODYSSEY vs VRT2")
    semilogy(alpha_vs_b[:,0],abs(vrt2_deflection_angle_int_fct(-alpha_vs_b[:,0])-alpha_vs_b[:,1]),'b+',label="ASTRORAY vs VRT2")
    # plot(impact_parameter_grid_retrograd,grtrans_deflection_angle_retrograd_int/pi,'k-',label="SPLINE")
    # plot(impact_parameter_grid_prograd,grtrans_deflection_angle_prograd_int/pi,'k-')
    axvspan(0-0.1,0+0.1,fc="k",alpha=0.5)
    axvspan(b_sc(a,-1),b_sc(a,1),fc="k",alpha=0.5)
    # text(0,0.1,"shadow",color="k",rotation="vertical",fontsize=20)
    text(0.51,0.51,"shadow",color="k",rotation="vertical",fontsize=20,transform=gca().transAxes,verticalalignment='center')
    grid(True)
    legend(labelspacing=0.1); setp(gca().get_legend().get_texts(), fontsize='12')
    axis((-20,20,None,None))
    xlabel(u"impact parameter "+r"$/M$")
    ylabel(r"$\Delta$"u" deflection angle / rad")


    # axis((None,None,0,1))

    tight_layout(h_pad=0.1,w_pad=0)

    savefig("deflection_angle_vs_impact_parameter.png")


if anim:
  figure(2)
  dangle=1
  for angle in range(angle_min, angle_max, dangle):
    ax3d.view_init(30, angle)
    magnify=cosh((angle-180.)/180.)/cosh(1)
    # magnify=magnify**10. # for camera-BH distance ~500rs this power gives good zoom range
    magnify=magnify**12. # for camera-BH distance ~500rs this power gives good zoom range
    # LIM=500*magnify
    # LIM=20 # 3
    # LIM=max(bmax,amax(r)*magnify) # RG:WIP Wanna take maximum b among all geodesics, not b from the last geodesic file in the list...
    LIM=max(10.,amax(r)*magnify) # RG:WIP Wanna take maximum b among all geodesics, not b from the last geodesic file in the list...

    ax3d.set_xlim3d(-LIM,LIM)
    ax3d.set_ylim3d(-LIM,LIM)
    ax3d.set_zlim3d(-LIM,LIM)
    if angle%100==0:
        print "drawing angle =",angle,"magnify =",magnify
    # draw()
    savefig("geodesics-3D-angle_"+string.zfill(angle,4)+".png")
    #savefig("geodesics-3D-angle_"+string.zfill(angle,4)+"_BADremoteGEODESICS.png")
    #savefig("geodesics-3D-angle_"+string.zfill(angle,4)+"_accur1e-2_step1e-2.png")
    #savefig("geodesics-3D-angle_"+string.zfill(angle,4)+"_accur1e-2.png")
    #savefig("geodesics-3D-angle_"+string.zfill(angle,4)+"_accur1e-6.png")

print "="*12+"\n====DONE====\n"+"="*12+"\n"
## DONE ##

