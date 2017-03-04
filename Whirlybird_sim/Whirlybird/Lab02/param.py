# Whirlybird Parameter File
import numpy as np

# From parameters list of the Whirlybird lab documentation
g = 9.81       # Gravity, m/s**2
L1 = 0.85      # Length of rod from vertex to the rotors, m
L2 = 0.3048    # Length of rod from vertex to counterbalance, m
m1 = 0.891     # Mass of weight at the end of l1 (between the rotors), kg
m2 = 1.00      # Mass of counterbalance at the end of l2, kg
d = 0.178      # Distance between m1 and a rotor (half the distance between rotors), m
h = 0.65       # Height of support rod, m
r = 0.12       # Length of one of the square rotors (maybe radius instead?), m
Jx = 0.0047    # Some force in x direction???, kg-m^2
Jy = 0.0014    # Some force in y direction???, kg-m^2
Jz = 0.0041    # Some force in z direction???, kg-m^2
#km = ????     # On parameters list, will add in later
#omega_gyro    # "
#omega_pixel   # "

# Other parameters for the animation
ell =  5      # For the axis, m

# Simulation Parameters
Ts = 0.01
sigma = 0.05

# For conversion to pwm
pwmConversionFactor = 10.875
pwm_e = .48
km = (m1*L1*g - m2*L2*g) / (L1 *(2*pwm_e))

# Initial Conditions
phi0 = 0.0*np.pi/180    # Pitch of Whirlybird relative to the ground, rads
phidot0 = 0.0           # Derivative of the pitch
theta0 = 0.0*np.pi/180  # Roll of Whirlybird relative to the ground, rads
thetadot0 = 0.0         # Derivative of the roll
psi0 = 0.0*np.pi/180    # Yaw of Whirlybird relative to the ground, rads
psidot0 = 0.0           # Derivative of the yaw

F_e = (m1 - m2*(L2/L1))*g
phimax = 45.0 *np.pi / 180.0
Fmax = 50.0
taumax = 50.0
thetamax = 60.0*np.pi / 180.0
####################################################
#                 PID Control:  Longitudinal
####################################################
# Flong_max = 10.85

# Open Loop
# kp*b0/(S**2 + kd*b*S + kp*b)
th_b0 = (L1/(m1*L1**2+m2*L2**2+Jy))
th_a1 = 0.0
th_a0 = 0.0

# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

th_tr = 1.4           # Rise time, s
th_zeta = 0.7      # Damping Coefficient
th_wn = 2.2/th_tr     # Natural frequency

# S**2 + alpha1*S + alpha0
th_alpha1 = 2.0*th_zeta*th_wn
th_alpha0 = th_wn**2

# Gains
# b0*kp/[S**2 + (a1 + b0*kd)S + (a0 + b0*kp)]
th_kp = (th_alpha0-th_a0)/th_b0
th_kd = (th_alpha1-th_a1)/th_b0

th_ki = 0.1
th_windup = 2

####################################################
#                 PID Control:  Lateral
####################################################
#Flat_max = ?

#---------------------------------------------------
#                    Inner Loop: Roll
#---------------------------------------------------
# Open Loop
# (kp/Jx)/(S**2 + (kp/Jx)*S + (kp/Jx))
phi_b0 = 1/Jx
# these are not used for some reason?
phi_a1 = 0.0
phi_a0  = 0.0

# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

phi_zeta = 0.7       # Damping Coefficient
phi_tr = 0.2          # Rise time, s
phi_wn = 2.2/phi_tr     # Natural frequency

# S**2 + alpha1*S + alpha0
phi_alpha1 = 2.0*phi_zeta*phi_wn
phi_alpha0 = phi_wn**2

# Gains
# b0*kp/[S**2 + (a1 + b0*kd)S + (a0 + b0*kp)]
phi_kp = (phi_alpha0)/phi_b0
phi_kd = (phi_alpha1)/phi_b0
phi_DC = 1   #DC gain

phi_ki = 0.0

#---------------------------------------------------
#                    Outer Loop: Yaw
#---------------------------------------------------
M = 10.0    # bandwidth separation factor
# Open Loop
# kp*b0/(S**2 + kd*b*S + kp*b)
psi_b0 = (L1*F_e/(m1*L1**2+m2*L2**2+Jz))
psi_a1 = 0.0
psi_a0 = 0.0

# Desired Closed Loop tuning parameters
# S**2 + 2*zeta*wn*S + wn**2

psi_tr = M*phi_tr  # Rise time, s
psi_zeta = 0.707      # Damping Coefficient
psi_wn = 2.2/psi_tr  # Natural frequency


# S**2 + alpha1*S + alpha0
psi_alpha1 = 2.0*psi_zeta*psi_wn
psi_alpha0 = psi_wn**2

# Gains
# b0*kp/[S**2 + (a1 + b0*kd*IL_DC)S + (a0 + b0*kp*IL_DC)]
psi_kp = (psi_alpha0-psi_a0)/(phi_DC*psi_b0)
psi_kd = (psi_alpha1-psi_a1)/(phi_DC*psi_b0)

psi_ki = 0.1
psi_windup = 0.15

print('km: ', km)
print('th_kp: ', th_kp)
print('th_kd: ', th_kd)
print('phi_kp: ',phi_kp)
print('phi_kd: ',phi_kd)
print('phi_DC:', phi_DC)
print('psi_kp: ',psi_kp)
print('psi_kd: ',psi_kd)
