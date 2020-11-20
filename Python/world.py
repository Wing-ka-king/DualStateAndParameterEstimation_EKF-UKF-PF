
import numpy as np
from math import cos, sin, tan, pi, sqrt, atan2, fabs
import matplotlib.pyplot as plt

import pickle
# from agent import Agent

def rand(thres):
	return (1 + thres*2*(np.random.rand() - 0.5))


class World:
	def __init__(self, width, height, animation_frequency=200):
		self.width = width
		self.height = height
		self.animation_frequency = animation_frequency
		self.start = None

		self.agents = []
		self.landmarks = []
		self.frequency = 100
		self.dt = 1./self.frequency

	def add_landmark(self, pos):
		# x, y, color
		self.landmarks.append(pos)

	# def create_new_agent(self, color="red", state=(0,0,0), name="name"):
	# 	a = Agent(self.dt, color=color, state=state, name=name)
	# 	a.world = self
	# 	self.agents.append(a)
	# 	return a

	def add_agent(self, agent):
		agent.world = self
		self.agents.append(agent)

	def plot(self):
		plt.pause(1./self.animation_frequency)
		plt.clf()

		# plot border
		plt.plot([0, self.width, self.width, 0, 0], [0, 0, self.height, self.height, 0])

		#plot landmarks

		for l in self.landmarks:
			x, y, color, RFID = l
			plt.plot([x], [y], color)

		# Plot ground truth pos
		for a in self.agents:
			a.plot()


		plt.axis([0, self.width, 0, self.height])
		plt.grid(True)

		plt.legend(loc='upper right')


	def load(self, name):
		color = {
			"left": 'bs',
			"right": "rs",
			"path": "gs"
		}
		with open(name, "rb") as handle:
			data = pickle.load(handle,encoding="latin1")

			RFID = 0

			self.track = data

			for d in data["left"]:
				self.add_landmark((d[0], d[1], color["left"], RFID))
				RFID += 1
			for d in data["right"]:
				self.add_landmark((d[0], d[1], color["right"], RFID))
				RFID += 1
			for d in data["path"]:
				RFID += 1
				self.add_landmark((d[0], d[1], color["path"], RFID))
			
			return data

	def create_control(self, cycle, dt):
		pwm = 5
		# remove first
		path = self.track["path"][1:, 0:]

		# Cyclic
		path = np.append(path, [path[0]], axis=0)
		# path = np.append(path, path[1:5], axis=0)

		control = []
		x, y = path[0][0], path[0][1]

		self.start = np.matrix([[x, y, 0, 0, 0]]).T

		agent.reset(self.start)
		x = self.start
		for i in range(1, len(path)):
			# if i > 2:
			# 	break
			w = 0
			c = 0
			while sqrt((path[i][1] - x[1,0])**2 + (path[i][0] - x[0,0])**2) > 5 or fabs(w) < pi/2:
				# if i == 18 and c > 20:
				# 	break
				c+=1

				# if (w < -pi ):
				# 	w = 2*pi + w
				# w = w % (2*pi)

				w = w*rand(0.0)

				w = min(pi/2, max(-pi/2, w))
				# if i==18:
				# 	print w, atan2(path[i][1] - y, path[i][0] - x), th, x, y, path[i]

				u = np.matrix([[v, w]]).T
				x, y = agent.step(dt, x, u)
				control.append(u)
				
				w = atan2(path[i][1] - x[1,0], path[i][0] - x[0,0]) - x[2,0]

				w = (w + pi)%(2*pi) - pi

		return control
		
		

