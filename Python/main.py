
import matplotlib.pyplot as plt
import numpy as np

import math
from math import cos, sin, tan, pi, atan2, sqrt, fabs

from world import World


from bicycle import Bicycle

import EKF


R_SIM = np.diag([0.01, 0.01, 0.1*pi/180, 0.01, 0, 0]);
Q_SIM = np.diag([0.05, 0.05, 0.01, 0.01]);


R_MODEL = np.diag([0.01, 0.01, 0.2*pi/180, 0.01, 0.0, 0.0]);
Q_MODEL = np.diag([0.05, 0.05, 0.015, 0.015])
# x,y,theta,w,a

EKF.R = R_MODEL
EKF.Q = Q_SIM

# state = [x, y, theta, v, 1/L, k/m].Transpose
# u = [pwm, steer].Transpose


W = World(100, 100)
W.load("track2.pkl")

N = 6000;
control = np.random.rand(2, N) - 0.5;


state = np.matrix([[10, 20, 0, 0, 1, 50]]).T


cycle_sim = Bicycle(state = state, color="green", name="Original");
cycle_sim.process_noise(R_SIM)
cycle_sim.observation_noise(Q_SIM)
W.add_agent(cycle_sim)


state = np.matrix([[10, 20, 0, 0, 0.1, .1]]).T

cycle_model = Bicycle(state = state, color="red", name="Model");
cycle_model.process_noise(R_MODEL)
cycle_model.observation_noise(Q_MODEL)
W.add_agent(cycle_model)

target = 1
path = W.track["path"][0:, 0:]
L = len(W.track["path"])
vmax = 0.5


for i in range(N):
	# u = control[:, i:i+1]
	# steer = atan2(path[target][1] - cycle_sim.state[1,0], path[target][0] - cycle_sim.state[0,0]) - cycle_sim.state[2,0]

	# steer = (steer + pi)%(2*pi) - pi;
	# if(fabs(steer) > pi/2) and sqrt((path[target][1] - cycle_sim.state[1,0])**2 + (path[target][0] - cycle_sim.state[0,0])**2) < 5:
	# 	steer = 0.;
	# 	target = (target + 1)%L;
	# steer = min(pi/2, max(-pi/2, steer))
	# u = np.matrix([[vmax - cycle_sim.state[3,0]],[ steer]])
	u = np.matrix([[2. + 0.*np.random.randn() - cycle_sim.state[3,0]],[ 0.5*sin(0.005*i) + 0.*np.random.randn()]])

	# Move actual system
	nX, nZ = cycle_sim.move(u)

	xEst, Qt = EKF.run(cycle_model, u, nZ);

	# xEst, nY = cycle_model.step(cycle_model.dt, cycle_model.state, u)
	cycle_model.set(xEst)

	if(i%W.frequency == 0):
		W.plot()





# velocity
def summary(name, i, type = "X"):
	fig1, ax = plt.subplots()


	ax.set_title(name)
	for a in W.agents:
		if a.name != "Original" and type == 'Y':
			continue
		data = a.historyX.A;
		if(type == 'Y'):
			data = a.historyY.A
		vel = data[i, 0:]
		ax.plot(vel, color=a.color, label=a.name)

	ax.legend(loc='upper right')

summary("V", 3)
summary("L_iv", 4)
summary("K/m", 5)
# summary("Measured w", 2, type="Y")
# summary("Measured acce", 3, type="Y")

plt.show()
