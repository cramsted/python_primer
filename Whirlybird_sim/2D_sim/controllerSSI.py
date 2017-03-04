import sys
import numpy as np 
import param as P

class controllerSSI:
  ''' This class inherits other controllers in order to organize multiple controllers.'''

  def __init__(self):
      # Instantiates the SSI_ctrl object
      self.theta_SSICtrl = thetaSSI_ctrl(P.Klat,P.kilat,P.z0,P.theta0,P.tau_max)
      self.h_SSICtrl = hSSI_ctrl(P.Klong,P.kilong,P.h0,P.F_max)
      # K is the closed loop SS gains
      # kr is the input gain
      # y0 is the initial position of the state

  # Clips u at the limit
  def saturate(self,limit,u):
     if abs(u) > limit:
          u = limit*np.sign(u)
     return u

  def getForces(self,y_r,y):
      # y_r is the referenced input
      # y is the current state  
      z_r = y_r[1]
      h_r = y_r[0]
      z = y[0]
      theta = y[2]
      h = y[1]
      
      # h_tilde = h - P.h0
      # h_r_tilde = h_r - P.h0 
      # z_tilde = z - P.z0
      # z_r_tilde = z_r - P.z0
      # theta_tilde = theta - P.theta0

      tau_unsat = self.theta_SSICtrl.SSI_loop(z_r,z,theta)
      tau_sat = self.saturate(P.tau_max,tau_unsat)
      Fe = P.F_e/np.cos(theta)

      F = self.h_SSICtrl.SSI_loop(h_r,h,theta) + Fe
      return [F, tau_sat]


class thetaSSI_ctrl:
  def __init__(self,K,ki,z0,theta0,limit):
      self.zdot = 0.0              # Difference term
      self.thetadot = 0.0          # Difference term
      self.integrator = 0.0
      self.z_d1 = z0               # Last z term
      self.theta_d1 = theta0       # Last theta term
      self.error_d1 = 0.0
      self.K = K                   # Closed loop SS gains
      self.ki = ki                 # Input gain
      self.limit = limit           # Maxiumum force
      self.input_distrubance = 0.1


  def SSI_loop(self,z_r,z,theta):
 
      error = z_r - z

      # Update Differentiator
      a1 = (2*P.sigma - P.Ts)/(2*P.sigma+P.Ts)
      a2 = 2/(2*P.sigma+P.Ts)
      
      self.zdot = a1*self.zdot + a2*(z-self.z_d1)
      self.thetadot = a1*self.thetadot + a2*(theta-self.theta_d1)

      self.integrator += (P.Ts/2.0)*(error+self.error_d1)


      self.z_d1 = z
      self.theta_d1 = theta 
      self.error_dl = error

      # Construct the state
      x = np.matrix([[z-z_r],
                     [theta],
                     [self.zdot],
                     [self.thetadot]])

      # Compute the state feedback controller
      tau_unsat = - self.K*x - self.ki*self.integrator #+ self.input_distrubance

      tau_sat = self.saturate(tau_unsat)

      if self.ki !=0:
        self.integrator += P.Ts/self.ki*(tau_sat-tau_unsat)

      return tau_sat.item(0)

  def saturate(self,u):
    if abs(u) > self.limit:
      u = self.limit*np.sign(u)
    return u 

class hSSI_ctrl:
  def __init__(self,K,ki,h0,limit):
      self.hdot = 0.0          # Difference term
      self.integrator = 0.0
      self.h_d1 = h0          # Delayed state
      self.error_d1 = 0.0
      self.limit = limit           # Max torque
      self.K = K                   # SS gains
      self.ki = ki                 # Input gain

  def SSI_loop(self,h_r,h,theta): 

      error = h_r - h

      # Update Differentiator
      a1 = (2*P.sigma - P.Ts)/(2*P.sigma+P.Ts)
      a2 = 2/(2*P.sigma+P.Ts)

      self.hdot = a1*self.hdot \
                          + a2*(h-self.h_d1)

      self.integrator += (P.Ts/2.0)*(error+self.error_d1)


      # Update theta_d1
      self.h_d1 = h
      self.error_d1 = error 

      # Construct the state
      x = np.matrix([[h],
                     [self.hdot]])

      # Compute the state feedback controller

      # Fe = P.F_e/np.cos(theta)
      F_unsat = -self.K*x - self.ki*self.integrator #+ Fe

      F_sat = self.saturate(F_unsat)

      if self.ki !=0:
        self.integrator += P.Ts/self.ki*(F_sat-F_unsat)

      return F_sat.item(0) 

  def saturate(self,u):
    if abs(u) > self.limit:
      u = self.limit*np.sign(u)
    return u 