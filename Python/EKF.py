
import matplotlib.pyplot as plt
import numpy as np

import math
from math import cos, sin, tan, pi, atan2, sqrt

N_STATE = 6

# state = [x, y, theta, v, 1/L, k/m].Transpose
# u = [pwm, steer].Transpose


P = np.diag(np.array([1,1,0.02,0.01,1.0,1.0]))

Qt = P

def compute_jacobian(dt, state, u):

	pwm = u[0,0]
	steer = u[1,0]

	v = state[3,0]
	theta = state[2,0]


	L_inv = state[4,0]
	K_m = state[5,0]


	G = np.matrix([	[ 1, 0, -dt*v*sin(theta), dt*cos(theta), 0, 0],
					[ 0, 1, dt*v*cos(theta), dt*sin(theta), 0, 0],
					[ 0, 0, 1, dt*L_inv*tan(steer), dt*v*tan(steer), 0],
					[ 0, 0, 0, 1, 0, dt*pwm],
					[ 0, 0, 0, 0, 1, 0],
					[ 0, 0, 0, 0, 0, 1]])

	H = np.matrix([ [1, 0, 0, 0, 0, 0],
					[0, 1, 0, 0, 0, 0],
					# [0, 0, 1, 0, 0, 0],
					[0, 0, 0, 0, v*tan(steer), 0],
					[0, 0, 0, 0, 0, pwm]])

	Jp = np.matrix([[v*tan(steer)*dt, 0],
					[0, pwm*dt]])

	return G, H




def run(cycle, u, z):
	global Qt

	xp, zp = cycle.step(cycle.dt, cycle.state, u);

	dz = z - zp


	G, H = compute_jacobian(cycle.dt, cycle.state, u)

	Qtp = np.matmul(G, np.matmul(Qt, G.T)) + R;


	S = np.matmul(H, np.matmul(Qtp, H.T)) + Q
	S_inv = np.linalg.pinv(S)
	K = np.matmul(Qtp, np.matmul(H.T, S_inv))

	xEst = xp + np.matmul(K, dz)
	Qt = np.matmul( (np.identity(N_STATE) - np.matmul(K, H)) , Qt);

	# xEst[4:6,:] = xEst[4:6,:]*0.5 + 0.5*cycle.state[4:6,:]


	return xEst, Qt






