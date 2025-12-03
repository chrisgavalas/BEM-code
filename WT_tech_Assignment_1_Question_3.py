import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import trapezoid

blade_data = np.loadtxt('bladedat.txt')
cylinder_data = np.loadtxt('cylinder.txt')
FFA_60 = np.loadtxt('FFA-W3-600.txt')
FFA_48 = np.loadtxt('FFA-W3-480.txt')
FFA_36 = np.loadtxt('FFA-W3-360.txt')
FFA_301 = np.loadtxt('FFA-W3-301.txt')
FFA_241 = np.loadtxt('FFA-W3-241.txt')
tc_data_airfoil = np.array([100.0, 60.0, 48.0, 36.0, 30.1, 24.1])

R = blade_data[len(blade_data)-1][0]

AOA_data = np.zeros(len(cylinder_data))
r_data = np.zeros(len(blade_data))
c_data = np.zeros(len(blade_data))
b_data = np.zeros(len(blade_data))
tc_data = np.zeros(len(blade_data))
C_l_100 = np.zeros(len(cylinder_data))
C_l_60 = np.zeros(len(cylinder_data))
C_l_48 = np.zeros(len(cylinder_data))
C_l_36 = np.zeros(len(cylinder_data))
C_l_301 = np.zeros(len(cylinder_data))
C_l_241 = np.zeros(len(cylinder_data))
C_d_100 = np.zeros(len(cylinder_data))
C_d_60 = np.zeros(len(cylinder_data))
C_d_48 = np.zeros(len(cylinder_data))
C_d_36 = np.zeros(len(cylinder_data))
C_d_301 = np.zeros(len(cylinder_data))
C_d_241 = np.zeros(len(cylinder_data))

for i in range(len(blade_data)):
    r_data[i] = blade_data[i][0]
    c_data[i] = blade_data[i][1]
    b_data[i] = blade_data[i][2]
    tc_data[i] = blade_data[i][3]

for i in range(len(cylinder_data)):
    AOA_data[i] = cylinder_data[i][0]*np.pi/180
    C_l_100[i] = cylinder_data[i][1]
    C_l_60[i] = FFA_60[i][1]
    C_l_48[i] = FFA_48[i][1]
    C_l_36[i] = FFA_36[i][1]
    C_l_301[i] = FFA_301[i][1]
    C_l_241[i] = FFA_241[i][1]
    C_d_100[i] = cylinder_data[i][2]
    C_d_60[i] = FFA_60[i][2]
    C_d_48[i] = FFA_48[i][2]
    C_d_36[i] = FFA_36[i][2]
    C_d_301[i] = FFA_301[i][2]
    C_d_241[i] = FFA_241[i][2]
    
C_l_data = np.array([C_l_100, C_l_60, C_l_48, C_l_36, C_l_301, C_l_241])
C_d_data = np.array([C_d_100, C_d_60, C_d_48, C_d_36, C_d_301, C_d_241])

interp_C_l = RegularGridInterpolator((tc_data_airfoil, AOA_data), C_l_data, method='cubic')
interp_C_d = RegularGridInterpolator((tc_data_airfoil, AOA_data), C_d_data, method='cubic')

r_new = np.linspace(blade_data[0][0],R,50)
c_new = np.interp(r_new,r_data,c_data)
b_new = np.interp(r_new,r_data,b_data)
tc_new = np.interp(r_new,r_data,tc_data)

tol = 1e-6
V_cut_out = 25
V_cut_in = 4
V_0_rated = 11.432343737333841
optimal_theta_p_C_P = -0.0119
omega_max = 0.978274
B = 3
rho = 1.225
dx = 1e-5
P_rated = 10.64*1e+6
V_0 = np.linspace(V_cut_in, V_cut_out, 10)
theta_P_opt = optimal_theta_p_C_P
omega = omega_max
Theta_p_stall = np.linspace(theta_P_opt, -20*np.pi/180, 10)
Theta_p_feather = np.linspace(theta_P_opt, 30*np.pi/180, 10)
Theta_p_s = np.zeros(len(V_0))
Theta_p_f = np.zeros(len(V_0))
for i_V, V in enumerate(V_0):
    if V < V_0_rated:
        Theta_p_s[i_V] = theta_P_opt
    else:
        for j_theta_p, theta_p in enumerate(Theta_p_stall): 
            p_n = np.zeros(len(r_new))
            p_t = np.zeros(len(r_new))
            for j in range(len(r_new)-1):
                chord = c_new[j]
                r = r_new[j]
                beta = b_new[j]*np.pi/180
                sigma = (chord*B)/(2*np.pi*r)
                a = 0
                a_prime = 0
                f = 0.1
                k = 0
                
                while k == 0:
        
                    phi = np.arctan(((1-a)*V_0[i_V])/((1+a_prime)*omega*r))
                    
                    aoa = phi - (theta_p + beta)
                    
                    interp_point = np.array([[tc_new[j], aoa]])
                    
                    C_l = interp_C_l(interp_point)[0]
                   
                    C_d = interp_C_d(interp_point)[0]
            
                    F = (2/np.pi)*np.arccos(np.exp(-(B*(R-r))/(2*r*np.sin(abs(phi)))))
        
                    C_n = C_l*np.cos(phi) + C_d*np.sin(phi)
        
                    C_t = C_l*np.sin(phi) - C_d*np.cos(phi)
        
                    C_T = (((1-a)**2)*C_n*sigma)/(np.sin(phi)**2)
        
                    a_old = a
        
                    a_prime_old = a_prime
        
                    if a <= 0.33:
                        
                        a_star = (sigma*C_n)*(1-a)/(4*F*np.sin(phi)**2)
                        
                    else:
                        
                        a_star = C_T/(4*F*(1-(1/4)*(5-3*a)*a))
    
                    a_star_prime = (sigma*C_t)*(1+a)/(4*F*np.sin(phi)*np.cos(phi))
        
                    a = f*a_star + (1-f)*a
        
                    a_prime = f*a_star_prime + (1-f)*a_prime
        
                    if abs(a_old - a) < 1e-6 and abs(a_prime_old - a_prime) < 1e-6:
                        
                        k = 1
                #print(C_l)
                aoa = phi - (theta_p + beta)
            
                V_rel_sq = (V_0[i_V]*(1-a))**2 + (omega*r*(1+a_prime))**2
                
                C_n = C_l*np.cos(phi) + C_d*np.sin(phi)
    
                C_t = C_l*np.sin(phi) - C_d*np.cos(phi)
    
                p_n[j] = (1/2)*rho*V_rel_sq*chord*C_n
    
                p_t[j] = (1/2)*rho*V_rel_sq*chord*C_t
                
                F = (2/np.pi)*np.arccos(np.exp(-(B*(R-r))/(2*r*np.sin(abs(phi)))))
                #print('for r =', r, 'and chord =', chord, 'and theta_p =', theta_p)
                #print('phi =', phi)
                #print('a =', a)
                #print('a_prime =', a_prime)
                #print('p_t =', p_t)
                #print('p_n =', p_n)
                #print('F =', F)
                #print('angle of attack is', aoa)
            #print('--------------------------------------------------------------------')
            #print('p_n:',p_n)
            #print('p_t:',p_t)
            
            integral = trapezoid(r_new*p_t, r_new, dx)
            integral_T = trapezoid(p_n, r_new, dx)
            P = omega*B*integral
            Thrust = B*integral_T
            #print("for V_0 =", V, "and theta_p =", theta_p)
            #print(P)
            if P > 0 and P - P_rated < -1e+3:
                Theta_p_s[i_V] = theta_p
                print("for V_0 =", V, "Max Power achieved at","theta_p_stalling =", theta_p, "and is P_max =", P)
                break

#------------------------------pitch to feather-----------------------------------------
for i_V, V in enumerate(V_0):
    if V < V_0_rated:
        Theta_p_f[i_V] = theta_P_opt
    else:
        for j_theta_p, theta_p in enumerate(Theta_p_feather): 
            p_n = np.zeros(len(r_new))
            p_t = np.zeros(len(r_new))
            for j in range(len(r_new)-1):
                chord = c_new[j]
                r = r_new[j]
                beta = b_new[j]*np.pi/180
                sigma = (chord*B)/(2*np.pi*r)
                a = 0
                a_prime = 0
                f = 0.1
                k = 0
                
                while k == 0:
        
                    phi = np.arctan(((1-a)*V_0[i_V])/((1+a_prime)*omega*r))
                    
                    aoa = phi - (theta_p + beta)
                    
                    interp_point = np.array([[tc_new[j], aoa]])
                    
                    C_l = interp_C_l(interp_point)[0]
                   
                    C_d = interp_C_d(interp_point)[0]
            
                    F = (2/np.pi)*np.arccos(np.exp(-(B*(R-r))/(2*r*np.sin(abs(phi)))))
        
                    C_n = C_l*np.cos(phi) + C_d*np.sin(phi)
        
                    C_t = C_l*np.sin(phi) - C_d*np.cos(phi)
        
                    C_T = (((1-a)**2)*C_n*sigma)/(np.sin(phi)**2)
        
                    a_old = a
        
                    a_prime_old = a_prime
        
                    if a <= 0.33:
                        
                        a_star = (sigma*C_n)*(1-a)/(4*F*np.sin(phi)**2)
                        
                    else:
                        
                        a_star = C_T/(4*F*(1-(1/4)*(5-3*a)*a))
    
                    a_star_prime = (sigma*C_t)*(1+a)/(4*F*np.sin(phi)*np.cos(phi))
        
                    a = f*a_star + (1-f)*a
        
                    a_prime = f*a_star_prime + (1-f)*a_prime
        
                    if abs(a_old - a) < 1e-6 and abs(a_prime_old - a_prime) < 1e-6:
                        
                        k = 1
                #print(C_l)
                aoa = phi - (theta_p + beta)
            
                V_rel_sq = (V_0[i_V]*(1-a))**2 + (omega*r*(1+a_prime))**2
                
                C_n = C_l*np.cos(phi) + C_d*np.sin(phi)
    
                C_t = C_l*np.sin(phi) - C_d*np.cos(phi)
    
                p_n[j] = (1/2)*rho*V_rel_sq*chord*C_n
    
                p_t[j] = (1/2)*rho*V_rel_sq*chord*C_t
                
                F = (2/np.pi)*np.arccos(np.exp(-(B*(R-r))/(2*r*np.sin(abs(phi)))))
                #print('for r =', r, 'and chord =', chord, 'and theta_p =', theta_p)
                #print('phi =', phi)
                #print('a =', a)
                #print('a_prime =', a_prime)
                #print('p_t =', p_t)
                #print('p_n =', p_n)
                #print('F =', F)
                #print('angle of attack is', aoa)
            #print('--------------------------------------------------------------------')
            #print('p_n:',p_n)
            #print('p_t:',p_t)
            
            integral = trapezoid(r_new*p_t, r_new, dx)
            integral_T = trapezoid(p_n, r_new, dx)
            P = omega*B*integral
            Thrust = B*integral_T
            #print("for V_0 =", V, "and theta_p =", theta_p)
            #print(P)
            if P > 0 and P - P_rated < -1e+3:
                Theta_p_f[i_V] = theta_p
                print("for V_0 =", V, "Max Power achieved at","theta_p_feathering =", theta_p, "and is P_max =", P)
                break
#----------------------------plotting θ_p vs V_0-------------------------------------------------
plt.plot(V_0, (Theta_p_f*180)/np.pi)
plt.plot(V_0, (Theta_p_s*180)/np.pi)
plt.xlabel('V_0')
plt.ylabel('theta_p')
plt.grid()
plt.show()