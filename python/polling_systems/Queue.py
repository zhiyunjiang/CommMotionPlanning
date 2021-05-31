import numpy as np

class Queue:

	def __init__(self):
		self.waiting = []
		self.wait_times = []
		self.is_being_serviced = False
		self.service_start_times = []

	def add(self, time):
		self.waiting.append(time)

	def start_service(self, time):
		if len(self.waiting) >= 0:
			in_time = self.waiting.pop(0)
			self.wait_times.append(time - in_time)
			self.service_start_times.append(in_time)			
			self.is_being_serviced = True

	def service_complete(self):
		self.is_being_serviced = False

	def avg_wait(self):
		if len(self.wait_times) > 0:
			return np.mean(self.wait_times)
		else:
			return 0




			

