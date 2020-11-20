import numpy as np
from math import cos, sin, tan, pi, sqrt, atan2, fabs
import matplotlib as mpl
import matplotlib.pyplot as plt



class Bicycle():
	"""docstring for Bicycle"""
	def __init__(self, state, color = "black", name = "Unnamed"):
		self.color = color;
		self.name = name;

		self.dt = 0.01
		self.state = state;
		self.historyX = state;

		self.nx = len(state.A[0:, 0])

		self.R = np.diag(np.zeros(self.nx))


		self.C = np.matrix([[1,0,0,0,0,0],
							[0,1,0,0,0,0],
							# [0,0,1,0,0,0],
							[0,0,0,0,0,0],
							[0,0,0,0,0,0]]);

		self.D = np.matrix([[0,0,0,0,0,0],
							[0,0,0,0,0,0],
							# [0,0,0,0,0,0],
							[0,0,1,0,0,0],
							[0,0,0,1,0,0]]);
		self.ny = len(self.D[:, 0])


		self.historyY = np.matrix(np.zeros((self.ny, 1)));



	def step(self, dt, state, u):
		# state = [x, y, theta, v, 1/L, k/m].Transpose
		# u = [pwm, steer].Transpose
		pwm = u[0,0];
		steer = u[1,0]*0.1;

		theta = state[2,0];
		v = state[3,0];
		L_inv = state[4,0];
		k_by_m = state[5,0];

		h_T = np.matrix([[v*cos(theta), v*sin(theta), v*L_inv*tan(steer), k_by_m*pwm, 0, 0]]).T;

		noise = np.random.randn(self.nx, 1)
		noise = np.matmul(self.R, noise)

		nX = state + h_T*dt + noise;

		nX[2,0] = (nX[2,0] + pi)%(2*pi) - pi


		noise = np.random.randn(self.ny, 1)
		noise = np.matmul(self.Q, noise)
		Y = np.matmul(self.C, nX) + np.matmul(self.D, h_T) + noise;

		return nX,  Y;



	def process_noise(self, R):
		self.R = R

	def observation_noise(self, Q):
		self.Q = Q

	def reset(self, state):
		self.state = state
		self.historyX = state
		self.historyY = np.matrix(np.zeros((self.ny, 1)));

	def set(self, state):
		self.state = state
		self.historyX = np.hstack((self.historyX, state))

	def move(self, u):
		nX, Y = self.step(self.dt, self.state, u);
		self.set(nX);

		self.historyY = np.hstack((self.historyY, Y));

		return nX, Y

	def plot(self):
		l = 5.
		x, y, th = self.state[0,0], self.state[1,0], self.state[2,0];
		# plt.plot([x], [y], "ro", markersize=8)

		t = mpl.markers.MarkerStyle(marker=">")
		t._transform = t.get_transform().rotate_deg(th*180/pi)

		plt.scatter((x), (y), marker=t, s= 100, color=self.color)
		# plt.plot(x, y, marker=(3, 1, th*180/pi - 90), markersize=20, linestyle='None')
		# plt.plot([x, x - l*cos(th)], [y, y - l*sin(th)], color="red", linewidth=2)

		# print self.historyX
		historyX = self.historyX.A
		plt.plot(historyX[0, 0:], historyX[1, 0:], color=self.color, label=self.name)

		# for l in self.observed_landmark:
		# 	plt.plot([x, l[1][0]], [y, l[1][1]], color="black")


# state = np.matrix([[10, 20, 0, 1, 1.0, 5.0]]).T
# u = np.matrix([[1, 9]]).T


# cycle = Bicycle()

# cycle.reset(state)
# cycle.set(state+state)
# cycle.set(state*3)

# print cycle.historyX