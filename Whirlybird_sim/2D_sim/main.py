import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import param as P
from step import step_function
from signal_generator import Signals
from sim_plot import plotGenerator
from controllerSSI import controllerSSI

# The Animation.py file is kept in the parent directory,
# so the parent directory path needs to be added.
from dynamics import planarVTOLDynamics
from animation import planarVTOLAnimation

# Converts force and torque into the left and
# right forces produced by the propellers.
def convertForces(u):
	F = u[0]         # Force, N
	tau = u[1]       # Torque, Nm
	# Convert Force and Torque to fl and fr
	# fl is the force created by the left propeller
	# fr is the force created by the right propeller
	fl = 1.0/(2.0)*(F-tau/P.d)
	fr = 1.0/(2.0)*(F+tau/P.d)
	u = [fl,fr]
	# print("F: ", F)
	# print("tau: ", tau)
	# print("u: ", u)
	return u

t_start = 0.0   # Start time of simulation
t_end = 20.0    # End time of simulation
t_Ts = P.Ts     # Simulation time step
t_elapse = 0.7  # Simulation time elapsed between each iteration
t_pause = 0.01  # Pause between each iteration


sig_gen = Signals()                   # Instantiate Signals class
plotGen = plotGenerator()             # Instantiate plotGenerator class
simAnimation = planarVTOLAnimation()  # Instantiate Animate class
dynam = planarVTOLDynamics()          # Instantiate Dynamics class
ctrl = controllerSSI()				  # Instatiate PID controller class

t = t_start               # Declare time variable to keep track of simulation time elapsed

while t < t_end:

    # Get referenced inputs from signal generators
	ref_input = sig_gen.getRefInputs(t)
	ref_input = [step_function(t,.5,5),step_function(t,.5,2.5)]

	# The dynamics of the model will be propagated in time by t_elapse
	# at intervals of t_Ts.
	t_temp = t +t_elapse
	while t < t_temp:

		states = dynam.States()              # Get current states
		u = ctrl.getForces(ref_input,states) # Calculate the forces
		u_convert = convertForces(u)         # Convert the forces from F,T to fl and fr
		dynam.propagateDynamics(u_convert)   # Propagate the dynamics of the model in time
		t = round(t +t_Ts,2)                 # Update time elapsed

	plt.figure(simAnimation.fig.number) # Switch current figure to animation figure
	simAnimation.drawPlanarVTOL(        # Update animation with current user input
		dynam.Outputs())
	plt.pause(0.0001)


	# Organizes the new data to be passed to plotGen
	new_data = [[ref_input[0],states[1]],
				[ref_input[1],states[0]],
				[0,states[2]],
				u]

	plotGen.updateDataHistory(t, new_data)

	plt.figure(plotGen.fig.number)		# Switch current figure to plotGen figure
	plotGen.update_plots()              # Update the plot
	plt.pause(0.0001)

	# time.sleep(t_pause)


plt.figure(plotGen.fig.number)
plotGen.update_plots()
plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('done')
input()
