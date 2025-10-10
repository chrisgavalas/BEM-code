import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import trapezoid
#-------------------------input data------------------------------------------
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

r_new = np.linspace(blade_data[0][0],R,100)
c_new = np.interp(r_new,r_data,c_data)
b_new = np.interp(r_new,r_data,b_data)
tc_new = np.interp(r_new,r_data,tc_data)

B = 3
rho = 1.225
V_0 = 10.0
dx = 1e-5
Theta_p = np.linspace(-4*np.pi/180, 3*np.pi/180, 20)
Lamda = np.linspace(5, 10, 20)
C_P_mat = np.zeros((len(Lamda),len(Theta_p)))
C_T_mat = np.zeros((len(Lamda),len(Theta_p)))
#-----------------------BEM algorithm-----------------------------------------
for i_lamda, lamda in enumerate(Lamda):
    for j_theta_p, theta_p in enumerate(Theta_p): 
        omega = (lamda*V_0)/R
        p_n = np.zeros(len(r_new))
        p_t = np.zeros(len(r_new))
        for i in range(len(r_new)-1):
            chord = c_new[i]
            r = r_new[i]
            beta = b_new[i]*np.pi/180
            sigma = (chord*B)/(2*np.pi*r)
            a = 0
            a_prime = 0
            f = 0.1
            k = 0
            
            while k == 0:
    
                phi = np.arctan(((1-a)*V_0)/((1+a_prime)*omega*r))
                
                aoa = phi - (theta_p + beta)
                
                interp_point = np.array([[tc_new[i], aoa]])
                
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
    
                if abs(a_old - a) < 1e-5 and abs(a_prime_old - a_prime) < 1e-5:
                    
                    k = 1
            #print(C_l)
            aoa = phi - (theta_p + beta)
        
            V_rel_sq = (V_0*(1-a))**2 + (omega*r*(1+a_prime))**2
            
            C_n = C_l*np.cos(phi) + C_d*np.sin(phi)

            C_t = C_l*np.sin(phi) - C_d*np.cos(phi)

            p_n[i] = (1/2)*rho*V_rel_sq*chord*C_n

            p_t[i] = (1/2)*rho*V_rel_sq*chord*C_t
            
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
        #print(P)
        C_P = P/(0.5*rho*(V_0**3)*np.pi*R**2)
        C_T = Thrust/(0.5*rho*(V_0**2)*np.pi*R**2)
        C_P_mat[i_lamda][j_theta_p] = C_P
        C_T_mat[i_lamda][j_theta_p] = C_T
        #plt.plot(r_new,p_t)
        #print(C_P)
max_index_C_P = np.unravel_index(np.argmax(C_P_mat), C_P_mat.shape)
max_index_C_T = np.unravel_index(np.argmax(C_T_mat), C_T_mat.shape)
max_lamda_index_C_P, max_theta_p_index_C_P = max_index_C_P
max_lamda_index_C_T, max_theta_p_index_C_T = max_index_C_T

# Get the actual parameter values
optimal_lamda_C_P = Lamda[max_lamda_index_C_P]
optimal_theta_p_C_P = Theta_p[max_theta_p_index_C_P]
optimal_lamda_C_T = Lamda[max_lamda_index_C_T]
optimal_theta_p_C_T = Theta_p[max_theta_p_index_C_T]
#C_T_max = np.max(C_T_mat)
#max_C_P = C_P_mat[max_lamda_index, max_theta_p_index]

#print(f"Maximum C_P: {max_C_P:.4f}")
print(f"Optimal Lamda (tip speed ratio): {optimal_lamda_C_P:.2f}")
print(f"Optimal Theta_p (pitch angle): {optimal_theta_p_C_P:.4f} rad ({optimal_theta_p_C_P*180/np.pi:.2f}°)")
print(np.max(C_P_mat))

#print(f"Optimal Lamda (tip speed ratio): {optimal_lamda_C_T:.2f}")
#print(f"Optimal Theta_p (pitch angle): {optimal_theta_p_C_T:.4f} rad ({optimal_theta_p_C_T*180/np.pi:.2f}°)")
#print(C_T_max)
# Create a contour plot of C_P vs Lamda and Theta_p
Theta_p_deg = Theta_p * 180/np.pi  # Convert to degrees for plotting

plt.figure(figsize=(10, 6))
contour_C_P = plt.contourf(Theta_p_deg, Lamda, C_P_mat, 20, cmap='viridis')
plt.colorbar(contour_C_P, label='Power Coefficient (C_P)')
plt.plot(optimal_theta_p_C_P*180/np.pi, optimal_lamda_C_P, 'ro', markersize=10, label='Optimal Point')
plt.xlabel('Pitch Angle (degrees)')
plt.ylabel('Tip Speed Ratio (λ)')
plt.title('Power Coefficient Contour Plot')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

plt.figure(figsize=(10, 6))
contour_C_T = plt.contourf(Theta_p_deg, Lamda, C_T_mat, 20, cmap='viridis')
plt.colorbar(contour_C_T, label='Thrust Coefficient (C_T)')
#plt.plot(optimal_theta_p_C_T*180/np.pi, optimal_lamda_C_T, 'ro', markersize=10, label='Optimal Point')
plt.xlabel('Pitch Angle (degrees)')
plt.ylabel('Tip Speed Ratio (λ)')
plt.title('Thrust Coefficient Contour Plot')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
#----------------printing results---------------------------------------------
            #print('for r =', r, 'and chord =', chord, 'and theta_p =', theta_p)
            #print('phi =', phi)
            #print('a =', a)
            #print('a_prime =', a_prime)
            #print('p_t =', p_t)
            #print('p_n =', p_n)
            #print('F =', F)
            #print('angle of attack is', aoa)
        #print('--------------------------------------------------------------------')

    
        
        
        
        
    
    