import numpy as np

n   = 24
off_1 = 1
off_2 = 2
#v_a = list(range(0,n))
#v_b = list(range(0,n))

v_a = np.ones(n)+off_1
v_b = np.ones(n)+off_2

v_anp = np.array(v_a)
v_bnp = np.array(v_b)

v_c = np.dot(v_anp, v_bnp)

print("size(va) = ", v_anp.size)
print("size(vb) = ", v_anp.size)
print("va       = ", v_anp)
print("vb       = ", v_bnp)
print("c        = ", v_c)
