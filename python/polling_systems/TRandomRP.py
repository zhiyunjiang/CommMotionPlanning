from MarkovianRP import RandomRP


def run_tests():
	test_constructor()

def test_constructor():
	assert _test_constructor_negative_pi()
	assert _test_constructor_sum_too_high()
	assert _test_constructor_OK()

def _test_constructor_negative_pi():
	success = False
	bad_pi = [-0.2, 1.2]
	try:
		rp = RandomRP(bad_pi1)
	except Exception as e:
		success = True

	return success


def _test_constructor_sum_too_high():
	success = False
	bad_pi = np.array([0.2, 1.1])
	try:
		rp = RandomRP(bad_pi)
	except Exception as e:
		success = True

	return success

def _test_constructor_OK():
	success = True
	pi = [0.9, 0.1]
	try:
		rp = RandomRP(pi)
	except Exception as e:
		success = False

	return success




if __name__== "__main___":
	run_tests()
