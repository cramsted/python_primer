import numpy as np
import param as P



class WhirlybirdDynamics:

    def __init__(self):

        # Initial state conditions
        self.state = np.matrix([[P.phi0],       # Phi initial orientation
                                [P.theta0],     # Theta initial orientation
                                [P.psi0],       # Psi initial orientation
                                [P.phidot0],    # Phidot initial velocity
                                [P.thetadot0],  # Thetadot initial velocity
                                [P.psidot0]])    # Psidot initial velocity

    def propagateDynamics(self,u):
        # P.Ts is the time step between function calls.
        # u contains the force and/or torque input(s).

        # RK4 integration
        k1 = self.Derivatives(self.state, u)
        k2 = self.Derivatives(self.state + P.Ts/2*k1, u)
        k3 = self.Derivatives(self.state + P.Ts/2*k2,u)
        k4 = self.Derivatives(self.state + P.Ts*k3, u)
        self.state += P.Ts/6 * (k1 + 2*k2 + 2*k3 + k4)


    # Return the derivatives of the continuous states
    def Derivatives(self,state,u):

        # States and forces
        phi = state.item(0)
        theta = state.item(1)
        psi = state.item(2)
        phidot = state.item(3)
        thetadot = state.item(4)
        psidot = state.item(5)
        fl = u[0]                   # make sure fl is from our main file function
        fr = u[1]
        # import pdb; pdb.set_trace()
        # ctheta and stheta are used multiple times. They are
        # precomputed and stored in another variable to increase
        # efficiency.
        cphi = np.cos(phi);
        sphi = np.sin(phi);
        ctheta = np.cos(theta);
        stheta = np.sin(theta);
        cpsi = np.cos(psi);
        spsi = np.sin(psi);


        # For convenience
        m1 = P.m1
        m2 = P.m2
        L1 = P.L1
        L2 = P.L2
        Jx = P.Jx
        Jy = P.Jy
        Jz = P.Jz

        # The equations of motion.
        M = np.matrix([[Jx,           0,                                                                   -Jx*stheta],
                       [0,              (m1*(L1**2))+(m2*(L2**2))+(Jy*(cphi**2))+(Jz*(sphi**2)), (Jy-Jz)*sphi*cphi*ctheta],
                       [-Jx*stheta,   (Jy-Jz)*sphi*cphi*ctheta,                                        ((m1*(L1**2))+(m2*(L2**2))+(Jy*(sphi**2))+(Jz*(cphi**2)))*(ctheta**2)+(Jx*(stheta**2))]])

        C = np.matrix([[-(thetadot**2)*(Jz-Jy)*sphi*cphi+(psidot**2)*(Jz-Jy)*sphi*cphi*(ctheta**2)-thetadot*psidot*ctheta*(Jx-(Jz-Jy)*((cphi**2)-(sphi**2)))],
                       [(psidot**2)*stheta*ctheta*(-Jx+m1*(L1**2)+m2*(L2**2)+Jy*(sphi**2)+Jz*(cphi**2))-2*phidot*thetadot*(Jz-Jy)*sphi*cphi-phidot*psidot*ctheta*(-Jx+(Jz-Jy)*((cphi**2)-(sphi**2)))],
                       [(thetadot**2)*(Jz-Jy)*sphi*cphi*stheta-phidot*thetadot*ctheta*(Jx+(Jz-Jy)*((cphi**2)-(sphi**2)))\
                        -2*phidot*psidot*(Jz-Jy)*(ctheta**2)*sphi*cphi+2*thetadot*psidot*stheta*ctheta*(Jx-m1*(L1**2)-m2*(L2**2)-Jy*(sphi**2)-Jz*(cphi**2))]])


        deltaP_deltaQ = np.matrix([[0],
                                   [(P.m1*P.L1-P.m2*P.L2)*P.g*ctheta],
                                   [0]])
        Q = np.matrix([[P.d*(fl-fr)],
                       [P.L1*(fl+fr)*cphi],
                       [P.L1*(fl+fr)*ctheta*sphi+P.d*(fr-fl)*stheta]])

        tmpM = np.linalg.inv(M)
        secondPart = Q-C-deltaP_deltaQ
        qddot = tmpM*secondPart

        phiddot = qddot.item(0)
        thetaddot = qddot.item(1)
        psiddot = qddot.item(2)

        xdot = np.matrix([[phidot],[thetadot],[psidot],[phiddot],[thetaddot],[psiddot]])

        return xdot


    # Returns the observable states
    def Outputs(self):
        # Return them in a list and not a matrix
        return self.state[0:6].T.tolist()[0]

    # Returns all current states
    def States(self):
        # Return them in a list and not a matrix
        return self.state.T.tolist()[0]
