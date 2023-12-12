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
A_p = 3
#f_S = 6E3
#A_S = 25
n = 1
B = 3E3
w0 = 10E3
#epsilon = sqrt(10**(A_p/10)-1)

w1 = 8611.87420807834 # in Hz
w2 = 11611.8742080783 # in Hz

# Print the values
print("w1 =", w1) # in Hz
print("w2 =", w2) # in Hz

w1_rads = w1*2*np.pi 
w2_rads = w2*2*np.pi
print("w1_rads =", w1_rads) # in rads
print("w2_rads =", w2_rads) # in rads

# Prototype Filter using the signal.cheby1

Wn = [w1_rads, w2_rads] 
rp = A_p
N=n
#num, den = signal.butter(n, f_p*2*np.pi/(epsilon**(1/n)),'low', analog=True)
num, den = signal.cheby1(N, rp , Wn, btype='bandpass', analog=True, output='ba')
w, h = signal.freqs(num, den, 1000) 

print('\n Den and num: \n')
print(den)
print(num)
print('\n\n')
# Transfer Function

R_1, R_2, R_3, C_1, C_2, R_A, R_B, R, C, alpha = symbols('R_1 R_2 R_3 C_1 C_2 R_A R_B R C alpha',positive = True)
v_in, v_x, v_out, v_y = symbols('v_in v_x v_out v_y')
s = sympy.Symbol('s')

# CHECK EQUATIONS
no1 = -(1/R_1 + 1/R_3 + s*C_1 + s*C_2)*v_x + (1/R_1)*v_in + (s*C_2)*v_y + (1/R_3)*v_out
no2 = (v_y-v_x)*s*C_2 + v_y/R_2
no3 = v_y*(1/R_A + 1/R_B) - v_out*(1/R_B) 

x = solve([no1,no2,no3],[v_x,v_y,v_out])

H = simplify(x[v_out]/v_in)
n,d = fraction(H)
simplify(d)
print('\n Simplified H: \n')
sympy.pretty_print(H)
print('\n\n')
#sympy.latex(H)

simplify(x[v_out]/x[v_y])

H1 = limit(H.subs([(R_A,R_A),(R_B,0)]),R_A,oo) #substitute R_A oo and R_B with 0
print('\n Simplified H1: \n')
sympy.pretty_print(H1)
print('\n\n')

# Define values
H2 = simplify(H1.subs([(R_1,R),(R_2,R),(R_3,R*alpha),(C_1,C),(C_2,C)]))

print('\n Simplified H2 final: \n')
sympy.pretty_print(H2)
print('\n\n')

# Find circuit values

n,d = fraction(H2)
simplify(d)
dd = poly(d,s)
aa = dd.coeffs()

p = roots(d,s)
print(f"\n roots: {list(p.values())}")
pp = list(p.keys())


#x = solve(#to do) # sub values
R1 = x[0][1]
R2 = x[0][1]
R3 = x[0][1]
C1 = x[0][0]
C2 = x[0][0]*x[0][2]
print('\n R1 = ', R1,'R2 = ', R2,'R3 = ', R3,'C1 = ',C1,'C2 = ', C2)



Hln = lambdify((s,R_1,R_2,R_3,C_1,C_2),H1)
Hln(sqrt(-1),R1,R2,R3,C1,C2)

x = np.abs(Hln(w*sqrt(-1),R1,R2,R3,C1,C2))
xx = np.asarray(x).astype(np.float64)
HldB = 20*np.log10(xx)

# Graph

plt.title('Cheby BP Sallen and Key filter frequency response')
plt.semilogx(w/(2*np.pi),20*np.log10(abs(h)), label = 'Prototype Cheby1 filter', color = 'dodgerblue')
plt.semilogx(w/(2*np.pi),HldB+0.2, label = 'Filter obtained from the circuit', color = 'peru')#+0.2??
plt.title('Butterworth prototype filter frequency response')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(w1, color ='green') # cutoff frequency
plt.axvline(w2, color ='red') # cutoff frequency
plt.legend()
plt.show()


