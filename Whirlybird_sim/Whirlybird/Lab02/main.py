import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import param as P
# from slider_input import Sliders
from signal_generator import Signals
from sim_plot import plotGenerator
from controlPD import controllerPID

# The Animation.py file is kept in the parent directory,
# so the parent directory path needs to be added.
sys.path.append('..')
from dynamics import WhirlybirdDynamics
from Animation import WhirlybirdAnimation

def convertForces(u):
		fl = .5*u[0]+(1/(2*P.d))*u[1]	   # This calculates fl from f and tau (from the sliders)
		fr = .5*u[0]-(1/(2*P.d))*u[1]	   # This calculates fr from f and tau (from the sliders)
		# import pdb; pdb.set_trace()
		propagateDerivativeIn = [fl,fr]	 # Store the forces in a list
		return propagateDerivativeIn		# The list will be passed in to propogateDerivative

def convertForcesToPWM(u):
	F = u[0]
	tau = u[1]
	pwml = (1/(2*P.km))*(F + tau / P.d)
	pwmr = (1/(2*P.km))*(F - tau / P.d)

	if pwml >= 0.6:
		pwml = 0.6
	if pwmr >= 0.6:
		pwmr = 0.6
	return [pwml,pwmr]

t_start = 0.0   # Start time of simulation
t_end = 15.0	# End time of simulation
t_Ts = P.Ts	 # Simulation time step
t_elapse = 0.1 # Simulation time elapsed between each iteration
t_pause = 0.01  # Pause between each iteration

plotGen = plotGenerator()
sig_gen = Signals()                   # Instantiate Signals class
simAnimation = WhirlybirdAnimation()  # Instantiate Animate class
dynam = WhirlybirdDynamics()			# Instantiate Dynamics class
ctrl = controllerPID()

t = t_start			   # Declare time variable to keep track of simulation time elapsed

while t < t_end:
    # Get referenced inputs from signal generators
	ref_input = sig_gen.getRefInputs(t)

	# The dynamics of the model will be propagated in time by t_elapse
	# at intervals of t_Ts.
	t_temp = t +t_elapse
	# import pdb; pdb.set_trace()

	while t < t_temp:

		states = dynam.Outputs()             # Get current states
		u = ctrl.getForces(ref_input,states) # Calculate the forces
		u_converted = convertForces(u)
		temp = convertForcesToPWM(u)
		# print("left: %.2f right %.2f" % (temp[0], temp[1]))
		dynam.propagateDynamics(u_converted)           # Propagate the dynamics of the model in time
		t = round(t +t_Ts,2)                 # Update time elapsed

	plt.figure(simAnimation.fig.number) # Switch current figure to animation figure
	simAnimation.drawWhirlybird(          # Update animation with current user input
		dynam.Outputs())
	plt.pause(0.0001)

	# Organizes the new data to be passed to plotGen
	new_data = [[ref_input[0], states[1]],
				[ref_input[0], states[0]],
				[ref_input[0], states[2]],
			    [temp[0],temp[1]]]
	plotGen.updateDataHistory(t, new_data)

	plt.figure(plotGen.fig.number)
	plotGen.update_plots()
	plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('done')
input()
# This function takes in the output u from the sliders F and tau
# and convert them into fl and fr as a list
