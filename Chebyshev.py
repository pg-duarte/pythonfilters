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
n = 2
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

Wn = [w1_rads, w2_rads] 
rp = A_p
N=n
#num, den = signal.butter(n, f_p*2*np.pi/(epsilon**(1/n)),'low', analog=True)
num, den = signal.cheby1(N, rp , Wn, btype='bandpass', analog=True, output='ba')
w, h = signal.freqs(num, den, 1000) 

# Prototype filter graph

plt.semilogx(w/(2*np.pi), 20 * np.log10(abs(h)), label = 'cheby1 filter')
plt.title('Butterworth prototype filter frequency response')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(w1, color ='green') # cutoff frequency
plt.axvline(w2, color ='red') # cutoff frequency
plt.legend()
plt.show()
print('k=',abs(h[1]),'num=',num,' denm',den)


