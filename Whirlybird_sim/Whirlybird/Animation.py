from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('Lab01/')
import param as P

class WhirlybirdAnimation:

	def __init__(self):
		self.flagInit=True   # Used to indicate initialization
		self.fig = plt.figure()
		self.ax = Axes3D(self.fig) # Create a 3D axes in the figure

		# A list object that will contain the lists of vertices of the
		# squares that represent the whirlybird's rotors.
		self.verts = self.getWhirlybirdVertices()

		# A list that will contain handles to the Poly3DCollection
		# Objects so that they can be modified.
		self.PolyCollections = []


		# Set axis limits
		_axis_limit = 1
		self.ax.set_zlim3d([-_axis_limit,_axis_limit])
		self.ax.set_ylim3d([-_axis_limit,_axis_limit])
		self.ax.set_xlim3d([-_axis_limit,_axis_limit])

		# Set title and labels
		self.ax.set_title('Whirlybird')
		self.ax.set_xlabel('East')
		self.ax.set_ylabel('North')
		self.ax.set_zlabel('-Down')

		# Change viewing angle
		self.ax.view_init(self.ax.elev, self.ax.azim+90)

	def getWhirlybirdVertices(self):

                # ADD 2 sticks here as L1 and L2 and H
		vertPoly = np.matrix([[0, P.L1, 0],					   #Center between L1 & L2
							  [-P.d, P.L1, 0],				   # Left rotar
							  [-P.d,       (P.L1+P.r/2), 0],   # Vertex 1 (X,Y,Z)
			      	          [(-P.d-P.r), (P.L1+P.r/2), 0],   # Vertex 2 (X,Y,Z)
			       	          [(-P.d-P.r), (P.L1-P.r/2), 0],   # Vertex 3 (X,Y,Z)
			       	          [-P.d,       (P.L1-P.r/2), 0],   # Vertex 4
    						  [-P.d, P.L1, 0],
  							  [P.d, P.L1, 0],				   # Right rotar
							  [P.d,       (P.L1+P.r/2), 0],    # Vertex 5 (X,Y,Z)
  			       	          [(P.d+P.r), (P.L1+P.r/2), 0],    # Vertex 6 (X,Y,Z)
  		       	              [(P.d+P.r), (P.L1-P.r/2), 0],    # Vertex 7 (X,Y,Z)
  			       	          [P.d,       (P.L1-P.r/2), 0],	   # Vertex 8
							  [P.d, P.L1, 0],
							  [0, P.L1, 0],					   # L2
							  [0, -P.L2, 0],				   # L2
							  [0, P.L1, 0]])
    # Vertices for the pole

    # Vertices for the vertical pole (should not rotate with the other vertices



		# Return the vertices in a list. verts 5 was unpacked to make the
		# rest of the code more computationally efficient, but it will be
		# packed back up later for demonstration purposes.
		return [vertPoly]

	def drawWhirlybird(self,u):
		# Desired configuration
		phi = u[0]
		theta = u[1]
		psi = u[2]
		# import pdb; pdb.set_trace()
		X = [0, 0]
		Y = [0, 0]
		Z = [0, -P.h]

		# Rotates and transforms the Whirlybird.
		verts = []
		for i in range(len(self.verts)):
			vertsTemp = rotate(self.verts[i].T,phi,theta,psi).T
			vertsTemp = transformXYZtoNED(vertsTemp)
			verts.append(vertsTemp)

		verts = np.asarray(verts)

		if self.flagInit == True:

	  		# Initialize Poly3DCollection class for each set of vertices, and
	  		# create an object handle to each one.
			self.PolyCollections.append(Poly3DCollection([verts[0]],
			            facecolor = 'red',edgecolor = 'black', lw = 2))
			line, =self.ax.plot(X, Y, Z, lw=3, c = 'black')


			# Add each Poly3DCollection object to the axes.
			for i in range(len(self.PolyCollections)):
				self.ax.add_collection3d(self.PolyCollections[i])
			self.flagInit = False
		else:
			# Update the verts
			for i in range(len(self.PolyCollections)):
				self.PolyCollections[i].set_verts([verts[i]])


# This function transforms the image from XYZ frame to NED frame
def transformXYZtoNED(XYZ):
	R = np.matrix([[0,1,0],
			[1,0,0],
			[0,0,-1]])
	NED = XYZ*R
	return NED

def rotate(XYZ,phi,theta,psi):

    # Define rotation matrix
    R_roll = np.matrix([[1,           0,            0],
            			[0, np.cos(-theta), -np.sin(-theta)],
            			[0, np.sin(-theta), np.cos(-theta)]])

    R_pitch = np.matrix([[ np.cos(phi), 0, np.sin(phi)],
            			 [0,              1,             0],
            			 [-np.sin(phi), 0, np.cos(phi)]])

    R_yaw = np.matrix([[np.cos(psi),-np.sin(psi), 0],
            		   [np.sin(psi), np.cos(psi), 0],
            		   [0,                     0, 1]])
    R = R_roll*R_pitch*R_yaw;

    # rotate vertices
    XYZ = R*XYZ

    return XYZ


if __name__ == "__main__":
	# Desired configuration of Whirlybird
	phi = 0.0*np.pi/180
	theta = 0.0*np.pi/180
	psi = 0.0*np.pi/180


	# Instantiate class
	draw_Whirlybird = WhirlybirdAnimation()

	# Draw the Whirlybird with desired configuration
	draw_Whirlybird.drawWhirlybird([phi,theta,psi])
	plt.show()
