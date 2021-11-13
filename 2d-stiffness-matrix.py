import numpy as np 


def calculate_S(E1, E2, v12, G12):

	S11 = 1 / E1
	S12 = -v12 / E1
	S22 = 1 / E2
	S66 = 1 / G12

	return S11, S12, S22, S66


def calculate_sin_cos(a):

	return np.sin(a), np.cos(a)


def calculate_Q(S11, S12, S22, S66):

	Q11 = S22 / (S11*S22 - S12**2)
	Q12 = -S12 / (S11*S22 - S12**2)
	Q22 = S11 / (S11*S22 - S12**2)
	Q66 = 1 / S66

	return Q11, Q12, Q22, Q66


def calculate_Qbar(Q11, Q12, Q22, Q66, a):

	s, c = calculate_sin_cos(a)

	Q11b = Q11*c**4 + Q22*s**4 + 2*(Q12 + 2*Q66)*s**2*c**2
	Q12b = (Q11 + Q22 - 4*Q66)*s**2*c**2 + Q12*(c**4 + s**2)
	Q22b = Q22*c**4 + Q11*s**4 + 2*(Q12 + 2*Q66)*c**2*s**2
	Q16b = (Q11 - Q12 - 2*Q66)*c**3*s - (Q22 - Q12 - 2*Q66)*s**3*c
	Q26b = (Q11 - Q12 - 2*Q66)*s**3*c - (Q22 - Q12 - 2*Q66)*c**3*s
	Q66b = (Q11 + Q22 - 2*Q12 - 2*Q66)*s**2*c**2 + Q66*(s**4 + c**4)

	return Q11b, Q12b, Q22b, Q16b, Q26b, Q66b


def make_stiffness_matrix(Q11b, Q12b, Q22b, Q16b, Q26b, Q66b):

	return np.array([[Q11b, Q12b, Q16b],
					 [Q12b, Q22b, Q26b],
					 [Q16b, Q26b, Q66b]])


def calculate_A(Qs, h):

	if len(Qs) != (len(h) - 1):
		return "The lengths are messed up!"

	A = np.zeros((3, 3))

	for count, ele in enumerate(Qs):

		A += Qs[count] * (h[count + 1] - h[count])

	return A


def calculate_B(Qs, h):

	if len(Qs) != (len(h) - 1):
		return "The lengths are messed up!"

	B = np.zeros((3, 3))

	for count, ele in enumerate(Qs):

		B += 0.5 * Qs[count] * (h[count + 1]**2 - h[count]**2)

	return B 


def calculate_D(Qs, h):

	if len(Qs) != (len(h) - 1):
		return "The lengths are messed up!"

	D = np.zeros((3, 3))

	for count, ele in enumerate(Qs):

		D += 1/3 * Qs[count] * (h[count + 1]**3 - h[count]**3)

	return D 


def main():

	# Alter this function based on your application

	# GPa   = 10**9
	E1    = 181
	v12   = 0.28
	E2    = 10.3
	G12   = 7.17
	th    = [0, np.pi/6, -np.pi/4]
	h     = [-0.0075, -0.0025, 0.0025, 0.0075]

	S11, S12, S22, S66 = calculate_S(E1, E2, v12, G12)
	Q11, Q12, Q22, Q66 = calculate_Q(S11, S12, S22, S66)

	Q = np.zeros((len(th), 3, 3))

	for i, a in enumerate(th):

		Q11b, Q12b, Q22b, Q16b, Q26b, Q66b = calculate_Qbar(Q11, Q12, Q22, Q66, a)
		Q[i] = make_stiffness_matrix(Q11b, Q12b, Q22b, Q16b, Q26b, Q66b)
		print(f"\n Reduced Transformed Stiffness Matrix Qbar for angle {np.degrees(a)} is (in GPa units): \n {Q[i]}")

	A = calculate_A(Q, h)
	B = calculate_B(Q, h)
	D = calculate_D(Q, h)

	print(f"\n The stiffness matrix A is (in GPa): \n {A}")
	print(f"\n The stiffness matrix B is (in GPa): \n {B}")
	print(f"\n The stiffness matrix D is (in GPa): \n {D}")


if __name__ == "__main__":

	main()

