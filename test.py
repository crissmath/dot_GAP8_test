import numpy as np





v_a = list(range(0,4))
v_b = list(range(0,4))

v_anp = np.array(v_a)
v_bnp = np.array(v_b)

v_c = np.dot(v_anp, v_bnp)
print("c = ", v_c)
