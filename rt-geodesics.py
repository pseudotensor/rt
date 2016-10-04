## DIAGNOSTIC ROUTINE TO ANALYZE GEODESICS COMPUTED WITH ASTRORAY

# USAGE:
# python -i rt-geodesics.py [USER-FILENAME] 0 0
# Last two numbers select range (min,max) in viewing angles

import string,matplotlib

try:
    __IPYTHON__ ## ARE WE RUNNING FROM IPYTHON?
except:
    matplotlib.use("Agg") # produces png without X

import pylab
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

angle_min=int(sys.argv[-2])
angle_max=int(sys.argv[-1])
sys.argv = sys.argv[:-2] # discard two last int arg for viewing angle

anim=False
figure1=True
figure4=True
figure3=True

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

for FILE in sys.argv[1:]:
    filename=FILE # e.g. ~> 93.75th43fn2150geodesic.dat

    d=loadtxt(filename)

    a=float(filename.split("geodesics")[1].split("th")[0])/100. # spin from filename
    rh=1.+sqrt(1.-a*a);


    lambda_affine=d[:,8]

    #WIP REPLACE
    # theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    # z = np.linspace(-2, 2, 100)
    # r = z**2 + 1
    # x = r * np.sin(theta)
    # y = r * np.cos(theta)

    # ASTRORAY [geoint.cpp]: y= t, phi, rp ,mup, r, mu ,tpr,phpr
    # outputting ppy[currth].cooxx[p][geo_idx] to file
    # BUT! evalpointzro only uses cooxx!=coox
    # coox[1][stN]=y[4]; # r
    # coox[2][stN]=y[5]; # mu=costh
    # coox[3][stN]=y[1]; # phi
    # ppy[currth].cooxx[p][stN]=coox[p][stN];
    
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

    tangent_in  = array((diff(x)[ 0],diff(y)[ 0],diff(z)[ 0]))
    tangent_out = array((diff(x)[-1],diff(y)[-1],diff(z)[-1]))

    # A.B=|A||B|cos(angle)
    cos_angle=dot(tangent_in,tangent_out)/norm(tangent_in)/norm(tangent_out)
    print "Deflection angle ",arccos(cos_angle),"rad =",arccos(cos_angle)/2./pi*360.,"deg"

    # http://arxiv.org/pdf/gr-qc/9907034v1.pdf eqs.(20),(24)
    # http://arxiv.org/pdf/1405.2919.pdf see Figs. 2.2,3.5
    # https://arxiv.org/pdf/0907.5352v1.pdf eq(34) exact bending angle in equatorial Kerr

    ##############################################################################################
    # !!! ===>   https://arxiv.org/pdf/gr-qc/0611086v2.pdf eqs.(8),(14->23),(17->25)    <=== !!! #
    ##############################################################################################
    r_per=amin(r) # For Schwarzschild b=r_per/sqrt(1.-r_s/r_per) http://arxiv.org/abs/gr-qc/9907034
    b=r_per # FIXME: wrong (not too bad for weak field though...)
    b_crit = 3.*sqrt(3.) # "photon sphere":3 https://arxiv.org/pdf/gr-qc/0611086v2.pdf eq 
    b_prime = 1.-b_crit/b # https://arxiv.org/pdf/gr-qc/0611086v2.pdf eq before eq (20) III.D page 12 ?

    # elliptic integrals: http://arxiv.org/pdf/math/9409227v1.pdf
    deflection_angle_0611086v2  = 0. # Schwarzschild need elliptic integral of the 1st kind

    try:
        from sympy.mpmath import ellippi as Pi
        EXACT_BENDING_ANGLE=True
    except ImportError:
        EXACT_BENDING_ANGLE=False
        continue

    # https://arxiv.org/pdf/0907.5352v1.pdf
    b_s = sign(a)*b # BH spin 
    # What about b_s=0 ?  :-s
    P = r_0*(1.+a/b_s)/(1.-a/b_s) # eq (16) in https://arxiv.org/pdf/0907.5352v1.pdf 
    Q = (P-2.)*(P+6.) # eq (19) in https://arxiv.org/pdf/0907.5352v1.pdf

    omega_s = a/b_s # eq (36)
    h = 1./r_0 # eq (36)
    h_sc = (1.+omega_s)/(1.-omega_s) # eq (37)
    omega_0 = a*a # eq (36)
    r_0_over_Q = 1. / h_sc / sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (40)  
    k_squared = (sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))+6.*h/h_sc-1.)/2./sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) # eq (41)
    psi = sqrt( ( 1.- 2.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) /  ( 1.- 6.*h/h_sc - sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) )) # eq (42)
    Omega_plus  =  +((1.+sqrt(1.-omega_0))*(1.-omega_s)-omega_0/2.)/sqrt(1.-omega_0)/(1.+sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)
    Omega_minus  =  -((1.-sqrt(1.-omega_0))*(1.-omega_s)+omega_0/2.)/sqrt(1.-omega_0)/(1.-sqrt(1.-omega_0)-omega_0*h_sc/4.*( 1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)) ) ) # eq (43)
    n_plus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.+sqrt(1.-omega_0))) # eq (43)
    n_minus = (1.-6*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc)))/(1.-2.*h/h_sc-sqrt((1.-2.*h/h_sc)*(1.+6.*h/h_sc))-4./omega_0/h_sc*(1.-sqrt(1.-omega_0))) # eq (43)

    if EXACT_BENDING_ANGLE==True:
        deflection_angle_09075352v1 = -pi+4./(1.-omega_s)*sqrt(r_0_over_Q)*(Omega_plus*(Pi(n_plus,k)-Pi(n_plus,psi,k))+Omega_minus*(Pi(n_minus,k)-Pi(n_minus,psi,k)))   # eq (34) in https://arxiv.org/pdf/0907.5352v1.pdf

    deflection_angle_weakfield   = 4./b * ( 1. + 15./16.*pi/b )
    deflection_angle_strongfield = log( 3.482/(b-3.*sqrt(3.)))
    print "pericenter=",r_per
    print "impact parameter b=",b
    print "...assuming Schwarzschild..."
    print "deflection_angle_weakfield=",deflection_angle_weakfield,"rad =",deflection_angle_weakfield/2./pi*360.,"deg"
    print "deflection_angle_strongfield=",deflection_angle_strongfield,"rad =",deflection_angle_strongfield/2./pi*360.,"deg"

    veto += size(find(r==None))
    if False: # (amax(r)>500 and size(lambda_affine[r==amax(r)])>1):
        veto += 1
    else:
        ax3d.plot(x, y, z,"b-",alpha=0.25,linewidth=1)


    #############
    # figure(1) #
    #############
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
# def drawSphere(xCenter, yCenter, zCenter, r):
#     #draw sphere
#     u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#     x=np.cos(u)*np.sin(v)
#     y=np.sin(u)*np.sin(v)
#     z=np.cos(v)
#     # shift and scale sphere
#     x = r*x + xCenter
#     y = r*y + yCenter
#     z = r*z + zCenter
#     return (x,y,z)

(xs,ys,zs) = drawSphere(0.,0.,0.,rh)
#ax3d.plot_surface(xs,ys,zs,color='k',alpha=1)
ax3d.plot_surface(xs,ys,zs,color='k')
ax3d.set_aspect("equal")

#################
## ANNOTATIONS ##

LIM=3
ax3d.set_xlim3d(-LIM,LIM)
ax3d.set_ylim3d(-LIM,LIM)
ax3d.set_zlim3d(-LIM,LIM)

#ax3d.dist = 5 # distance to camera DEF:10 
#ax3d.legend()

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


if anim:
  figure(2)
  for angle in range(angle_min, angle_max, 1):
    ax3d.view_init(30, angle)
    magnify=cosh((angle-180.)/180.)/cosh(1)
    # magnify=magnify**10. # for camera-BH distance ~500rs this power gives good zoom range
    magnify=magnify**12. # for camera-BH distance ~500rs this power gives good zoom range
    LIM=500*magnify
    # LIM=3
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

