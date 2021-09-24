class RRange():

	def __init__(self, r_max, r_min):
		self.max = r_max
		self.min = r_min

	def __str__(self):
		return 'max: %f, min: %f'%(self.max, self.min)

	def includes(self, r):
		return self.min <= r <= self.max