#import bodas
import sympy
from sympy import *

import numpy as np
import math

from scipy import signal
import matplotlib.pyplot as plt

import control

#Low Pass Butterworth Filter Requirements:
f_p = 1E3
A_p = 1
#f_S = 6E3
#A_S = 25
n = 2
epsilon = sqrt(10**(A_p/10)-1)
"""
for n in range(1,10):
    if 10*np.log10(float(1+(epson**2)*(f_S/f_p)**(2*n)))>A_S:
        break
print('N =',n,'As =',10*np.log10(float(1+(epson**2)*(f_S/f_p)**(2*n))))
"""
num, den = signal.butter(n, f_p*2*np.pi/(epsilon**(1/n)),'low', analog=True)
w, h = signal.freqs(num, den, 1000)

# Prototype filter graph
"""
plt.semilogx(w/(2*np.pi), 20 * np.log10(abs(h)))
plt.title('Butterworth prototype filter frequency response')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(f_p, color ='green') # cutoff frequency
plt.show()
print('k=',abs(h[1]),'num=',num,' denm',den)
"""

#Transfer Func
R_1, R_2, C_1, C_2, R_A, R_B, R, C, alpha = symbols('R_1 R_2 C_1 C_2 R_A R_B R C alpha',positive = True)
v_in, v_x, v_out, v_y = symbols('v_in v_x v_out v_y')
s = sympy.Symbol('s')

no1 = (1/R_1 + 1/R_2 + s*C_1)*v_x - (1/R_1)*v_in - (1/R_2)*v_y - s*C_1*v_out
no2 = (1/R_2 + s*C_2)*v_y - (1/R_2)*v_x
no3 = v_y*(R_B/R_A) + v_y - v_out

x = solve([no1,no2,no3],[v_x,v_y,v_out])

print(x) # clean
H = simplify(x[v_out]/v_in)
n,d = fraction(H)
simplify(d)
#sympy.pretty_print(H)
#sympy.latex(H)

simplify(x[v_out]/x[v_y])

H1 = limit(H.subs([(R_A,R_A),(R_B,0)]),R_A,oo)
#sympy.pretty_print(H1)

# Define values
H2 = simplify(H1.subs([(R_1,R),(R_2,R),(C_1,C),(C_2,alpha*C)]))

print('\n Simplified: \n')
sympy.pretty_print(H2)
print('\n\n')

n,d = fraction(H2)
simplify(d)
dd = poly(d,s)
aa = dd.coeffs()
wp = simplify(sqrt(aa[2]/aa[0]))
Qp = simplify(wp*aa[0]/aa[1])

print('wp:')
sympy.pretty_print(wp)
print('\n Qp:')
sympy.pretty_print(Qp)

p = roots(d,s)
print(f"\n roots: {list(p.values())}")
pp = list(p.keys())

k = num[0]/den[2]

x = solve([den[1]-aa[1]/aa[0],R-10e3,den[2]-aa[2]/aa[0]],[C,R,alpha])
R1 = x[0][1]
R2 = x[0][1]
C1 = x[0][0]
C2 = x[0][0]*x[0][2]
print('\n R1 = ', R1,'R2 = ', R2,'C1 = ', C1,'C2 = ', C2,)

Hln = lambdify((s,R_1,R_2,C_1,C_2),H1)
Hln(sqrt(-1),R1,R2,C1,C2)

x = np.abs(Hln(w*sqrt(-1),R1,R2,C1,C2))
xx = np.asarray(x).astype(np.float64)
HldB = 20*np.log10(xx)

plt.semilogx(w/(2*np.pi),20*np.log10(abs(h)), label = 'Prototype filter', color = 'dodgerblue')
plt.semilogx(w/(2*np.pi),HldB+0.2, label = 'Filter obtained from the circuit', color = 'peru')#+0.2??

plt.title('Sallen and Key filter frequency response')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(f_p, color ='green') # cutoff frequency
plt.legend()
plt.show()

