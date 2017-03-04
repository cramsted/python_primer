import sys
import numpy as np
import param as P

class controllerPID:
  ''' This class inherits other controllers in order to organize multiple controllers.'''

  def __init__(self):
      # Instantiates the PD_ctrl object
      self.psiCtrl = psiPID_ctrl(P.psi_kp,P.psi_kd,P.psi_ki,P.psi0,P.phimax,
        P.psi_windup)
      self.phiCtrl = phiPID_ctrl(P.phi_kp,P.phi_kd,P.phi_ki,P.phi0,P.taumax)
      self.thetaCtrl = thetaPID_ctrl(P.th_kp,P.th_kd,P.th_ki,P.theta0,
        P.thetamax,P.th_windup)
      # kp is the proportional gain
      # kd is the derivative gain
      # ki is the integral gain
      # y0 is the initial position of the state
      # P.error_max is the maximum error before saturation

  # Clips u at the limit
  def saturate(self,limit,u):
     if abs(u) > limit:
          u = limit*np.sign(u)
     return u


  def getForces(self,y_r,y):
    # y_r is the referenced input
    # y is the current state
    theta = y[1]
    psi = y[2]
    phi = y[0]
    theta_r = y_r[0]
    psi_r = y_r[1]
    thetadot = y[4]
    psidot = y[5]
    phidot = y[3]

    # F_t = P.th_kp * (theta_r - theta) - P.th_kd * thetadot # Calculate the force output
    Fe = (((P.m1*P.L1-P.m2*P.L2) * P.g) / P.L1) * np.cos(theta)
    F = self.thetaCtrl.thetaPID_loop(theta_r,theta) + Fe

    phi_r = self.psiCtrl.psiPID_loop(psi_r,psi)
    tau = self.phiCtrl.phiPID_loop(phi_r,phi)
    # if F > P.F_max:
    #     F = P.F_max
    return [F, tau]

class phiPID_ctrl:
  def __init__(self,kp,kd,ki,phi0,limit):
      self.phidot = 0.0    # Difference term
      self.integrator = 0.0        # Integrator term
      self.phi_d1 = phi0       # Delayed y output
      self.error_d1 = 0.0          # Delayed error
      self.kp = kp                 # Proportional control gain
      self.kd = kd                 # Derivative control gain
      self.ki = ki                 # Integral controal gain
      self.limit = limit           # Max torque

  def phiPID_loop(self,phi_r,phi):
      # Compute the current error
      error = phi_r - phi

      # Update phidot
      a1 = (2*P.sigma - P.Ts)/(2*P.sigma+P.Ts)
      a2 = 2/(2*P.sigma+P.Ts)
      self.phidot = a1*self.phidot \
                          + a2*(phi - self.phi_d1)
      self.phi_d1 = phi # update the delay

      # Update Integrator
      if abs(self.phidot) <0.01:
        self.integrator += (P.Ts/2.0)*(error+self.error_d1)

      # Update error_d1
      self.error_d1 = error

      # PD Control to calculate T
      tau_unsat = self.kp*error - self.kd*self.phidot + \
                self.ki*self.integrator

      tau_sat = self.saturate(tau_unsat)

    #   if self.integrator > 5:
    #     self.integrator = 5
      return tau_sat


  def saturate(self,u):
    if abs(u) > self.limit:
      u = self.limit*np.sign(u)
    return u


class psiPID_ctrl:
  def __init__(self,kp,kd,ki,z0,limit,windup):
      self.differentiator = 0.0    # Difference term
      self.integrator = 0.0        # Integral term
      self.psi_d1 = z0               # Delayed y output
      self.error_d1 = 0.0          # Delayed error
      self.kp = kp                 # Proportional control gain
      self.kd = kd                 # Derivative control gain
      self.ki = ki                 # Integral control gain
      self.limit = limit           # Maximum phi
      self.windup = windup


  def psiPID_loop(self,psi_r,psi):
      # Compute the current error
      error = psi_r - psi
      
      # Update Differentiator
      a1 = (2*P.sigma - P.Ts)/(2*P.sigma+P.Ts)
      a2 = 2/(2*P.sigma+P.Ts)
      self.differentiator = a1*self.differentiator \
                          + a2*(psi - self.psi_d1)

      self.psi_d1 = psi

      # Update Integrator
      # Only integrate when moving slowly
      if abs(self.differentiator) < 0.05:
        self.integrator += (P.Ts/2.0)*(error+self.error_d1)

      # UPIDate error_d1
      self.error_d1 = error

      # PID Control to calculate T
      phi_r_unsat = self.kp*error - self.kd*self.differentiator + self.ki*self.integrator

      phi_r_sat = self.saturate(phi_r_unsat)

      # print("integrator: ", self.integrator)
      if abs(self.integrator) > self.windup:
        self.integrator = self.windup*np.sign(self.integrator)

      return phi_r_sat

  def saturate(self,u):
    if abs(u) > self.limit:
      u = self.limit*np.sign(u)
    return u

class thetaPID_ctrl:
  def __init__(self,kp,kd,ki,z0,limit,windup):
      self.differentiator = 0.0    # Difference term
      self.integrator = 0.0        # Integral term
      self.theta_d1 = z0               # Delayed y output
      self.error_d1 = 0.0          # Delayed error
      self.kp = kp                 # Proportional control gain
      self.kd = kd                 # Derivative control gain
      self.ki = ki                 # Integral control gain
      self.limit = limit           # Maximum F
      self.windup = windup

  def thetaPID_loop(self,theta_r,theta):
      # Compute the current error
      error = theta_r - theta

      # Update Differentiator
      a1 = (2*P.sigma - P.Ts)/(2*P.sigma+P.Ts)
      a2 = 2/(2*P.sigma+P.Ts)
      self.differentiator = a1*self.differentiator \
                          + a2*(theta -self.theta_d1)

      self.theta_d1 = theta
      # Update Integrator
      # Only update integrator if moving slow
      if abs(self.differentiator) <0.05:
        self.integrator += (P.Ts/2.0)*(error+self.error_d1)

      # Update error_d1
      self.error_d1 = error

      # PID Control to calculate T
      F_r_unsat = self.kp*error - self.kd*self.differentiator + self.ki*self.integrator

      F_r_sat = self.saturate(F_r_unsat)

      if abs(self.integrator) > self.windup:
        self.integrator = self.windup*np.sign(self.integrator)

      return F_r_sat

  def saturate(self,u):
    if abs(u) > self.limit:
      u = self.limit*np.sign(u)
    return u

