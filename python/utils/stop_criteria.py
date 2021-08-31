class StopCriteria():

	def __init__(self, solution_required, max_iterations=float('inf'), max_run_time=float('inf')):
		self.solutionRequired = solution_required
		self.maxIterations = max_iterations
		self.maxRunTime = max_run_time

	def stop(self, sol_found, iterations, time, conj = 'OR'):
		if not sol_found and self.solutionRequired:
			return False

		if conj == 'OR':
			return (iterations >= self.maxIterations) or (time>= self.maxRunTime)
		else:
			return (iterations >= self.maxIterations) and (time>= self.maxRunTime)

