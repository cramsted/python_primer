import numpy as np
import control as cnt

# Physical parameters of the inverted pendulum
mc = 1.0        # Center mass, kg
Jc = 0.0042     # Inertia, kg m**2
mr = 0.25       # Right motor mass, kg
ml = 0.25       # Left motor mass, kg
d = 0.3         # Rod length, m
mu = 0.1        # Drag coeff, kg/s
g = 9.81        # Gravity m/s**s
h_max = 5       # Maximum Altitude Distance, m
z_max = 5       # Minimum Latitude Distance, m

# Simulation Parameters
Ts = 0.01
sigma = 0.05

# Initial Conditions
z0 = 0.0                # Latitude, m
h0 = 0.0                # Altitude, m
theta0 = 0.0*np.pi/180  # ,rads
zdot0 = 0.0
hdot0 = 0.0
thetadot0 = 0.0

zv0 = 0.0               # Position of target, m
includeTarget = False   # Determines whether to render target or no.

F_e = (mc + 2*mr) * g
fmax = 10
tau_max = (fmax-F_e/2)/d

mixing = np.array([[1,1],
				   [d,-d]])
####################################################
#                 State Space: lat
####################################################
F_max = 5                    # Max Force, N
error_max = 1                # Max step size,m
theta_max = 30.0*np.pi/180.0 # Max theta 
# M = 5                        # Time scale separation between
					         # inner and outer loop

# State Space Equations
# xdot = A*x + B*u
# y = C*x

Alat = np.matrix([[0.0,0.0,               1.0,      0.0],
			   [0.0,0.0,               0.0,      1.0],
			   [0.0,-F_e/(mc+2*mr),  -mu/(mc+2*mr),  0.0],
			   [0.0,0.0, 0.0, 0.0]])

Blat = np.matrix([[0.0],
			   [0.0], 
			   [0.0],
			   [1.0/(Jc+2*mr*d**2)]])

Clat = np.matrix([[1.0,0.0,0.0,0.0],
			   [0.0,1.0,0.0,0.0]])

Crlat = Clat[0,:]

# Augmented Matrices
A1lat = np.concatenate((
	np.concatenate((Alat,np.zeros((4,1))),axis=1),
	np.concatenate((-Crlat,np.matrix([[0.0]])),axis=1)),axis = 0)



B1lat = np.concatenate((Blat,np.matrix([[0.0]])),axis = 0)
# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

# import pdb; pdb.set_trace()
th_tr = 2.22           # Rise time, s
th_zeta = 0.707       # Damping Coefficient
th_wn = 2.2/th_tr     # Natural frequency

# S**2 + alpha1*S + alpha0
th_alpha1 = 2.0*th_zeta*th_wn
th_alpha0 = th_wn**2 


# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

z_tr = 5.0      # Rise time, s
z_zeta = 0.707      # Damping Coefficient
z_wn = 2.2/z_tr     # Natural frequency
lat_integrator_pole = -10.0

# S**2 + alpha1*S + alpha0
z_alpha1 = 2.0*z_zeta*z_wn
z_alpha0 = z_wn**2 

# Desired Poles
lat_des_char_poly = np.convolve(
	np.convolve([1,z_alpha1,z_alpha0],[1,th_alpha1,th_alpha0]),
	np.poly(lat_integrator_pole))
lat_des_poles = np.roots(lat_des_char_poly)

# Controllability Matrix
if np.linalg.matrix_rank(cnt.ctrb(A1lat,B1lat))!=5:
	print("The system is not controllable")
else:
	K1lat = cnt.acker(A1lat,B1lat,lat_des_poles)
	Klat = K1lat[0,0:4]
	kilat = K1lat[0,4]
	print('K1lat: ', K1lat)
	print('Klat: ', Klat)
	print('kilat: ', kilat)
	# print('Crlat: ', Crlat)


####################################################
#                 State Space: long
###################################################

# State Space Equations
# xdot = A*x + B*u
# y = C*x

Along = np.matrix([[0.0,1.0],
			       [0.0,0.0]])

Blong = np.matrix([[0.0],
			   [1/(mc+2*mr)]])

Clong = np.matrix([[1.0,0.0]])

# Form augmented system
A1long = np.matrix([[0.0,1.0,0.0],
			    [0.0,0.0,0.0],
			    [-1.0,0.0,0.0]])
B1long = np.matrix([[0.0],
				[Blong.item(1)],
				[0.0]])

# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

h_tr = 2.2      # Rise time, s
h_zeta = 0.707      # Damping Coefficient
h_wn = 2.2/h_tr     # Natural frequency
long_integrator_pole = -3.0

# S**2 + alpha1*S + alpha0
h_alpha1 = 2.0*h_zeta*h_wn
h_alpha0 = h_wn**2 

# Desired Poles
long_des_char_poly = np.convolve([1,h_alpha1,h_alpha0],np.poly(long_integrator_pole))
long_des_poles = np.roots(long_des_char_poly)

# Controllability Matrix
if np.linalg.matrix_rank(cnt.ctrb(A1long,B1long))!=3:
	print("The system is not controllable")
else:
	K1long = cnt.acker(A1long,B1long,long_des_poles)
	Klong = K1long[0,0:2]
	kilong = K1long[0,2]
	print('K1long: ', K1long)
	print('Klong: ', Klong)
	print('kilong: ',kilong)


#################################################
#          Uncertainty Parameters
#################################################
UNCERTAINTY_PARAMETERS = False
if UNCERTAINTY_PARAMETERS:
    alpha = 0.2;                                  # Uncertainty parameter
    mc = mc*(1+2*alpha*np.random.rand()-alpha)  # Mass
    Jc = Jc*(1+2*alpha*np.random.rand()-alpha)   # spring constant
    d = d*(1+2*alpha*np.random.rand()-alpha)   # frictional force
    mu = mu*(1+2*alpha*np.random.rand()-alpha)   # frictional force

print('mc: ', mc)
print('Jc: ', Jc)
print('d: ', d)
print('mu: ', mu)