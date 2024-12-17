import mean_val as m
import matplotlib.pyplot as plt
import numpy as np
l_s = [213.1, 187.1, 179.4, 147.0, 113.7] #cm
n_s = np.array([1, 2,3,4,5])

#y = ax+b
# rn(l) = a*l + rn(0)

#y = rn(l), a = a,  x = l, b = rn(0)
#njegov nulti moj prvi
#r1
r_1l = [0.45, 0.4, 0.4, 0.35, 0.3]

#r2 
r_2l = [1.0, 0.9, 0.9, 0.8, 0.75]

#r3
r_3l = [1.4, 1.3, 1.2, 1.15, 1.1]

#r4
r_4l = [1.7, 1.52, 1.5, 1.4, 1.3]

#r5
r_5l = [1.9, 1.85, 1.8, 1.6, 1.52]

rl_s = np.array([r_1l, r_2l, r_3l, r_4l, r_5l])

koef_a = []
koef_b = []
dev_a =[]
dev_b = []

for el in n_s:
    mr = m.Value(l_s, rl_s[el-1])
    a,b = mr.mean_v()
    deva, devb = mr.devijacija_a_b()
    koef_a.append(a)
    koef_b.append(b)
    dev_a.append(deva)
    dev_b.append(devb)

r0_s = np.array(koef_b)**2
#y = rn^2 , x = n, A = R lambda , b = -2rd0
mj_2 = m.Value(n_s, r0_s)
k1,k2 = mj_2.mean_v()
g1,g2 = mj_2.devijacija_a_b()
#av = k1/10000 #cm^2 u m^2
#print( (av/632.8 )* 10**9)
#print(((g1*10**(-4))/632.8 )* 10**9)

print(k1, k2)
def f(x):
    a = k1
    b=k2
    return a*x + b

plt.scatter(n_s, r0_s, label= 'mjerenja',  color = 'red')
plt.plot(n_s, f(n_s), label= 'regresijski pravac', lw = 0.7)
#plt.title(f'Leća: {} mm')
plt.xlabel('n')
plt.ylabel(r'r$^{2}_{0}$ [cm$^{2}$]')
plt.legend(loc = 'upper left')
plt.grid(True) 

plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.show()