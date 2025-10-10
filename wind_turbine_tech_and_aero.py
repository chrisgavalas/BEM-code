import numpy as np
import scipy as sp
from scipy.optimize import fsolve
from Glauert import Glauert
R = 1 
B = 3 
aoa = 4*(np.pi/180)
C_l = 1.201
C_d = 0.014
#TSR = 0.7 #tipspeedratio
lamda = 8
F = 1
TSR = np.linspace(0.7, 1.0, 10)
for tsr in TSR:
    x = lamda*tsr
    result = Glauert(R, B, aoa, C_l, C_d, lamda, F, x)
    theta = result[0]
    a = result[1]
    local_chord = result[2]
    
    print('for TSR =', tsr, 'the results are')
    print ('theta is:', theta)
    print('a is:', a)
    print('local_chord is:', local_chord)
    
    
    
    




