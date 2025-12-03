import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import trapezoid
from deflections import deflections
from numpy.linalg import eig
#-------------------------input data------------------------------------------
blade_data = np.loadtxt('bladedat.txt')
cylinder_data = np.loadtxt('cylinder.txt')
FFA_60 = np.loadtxt('FFA-W3-600.txt')
FFA_48 = np.loadtxt('FFA-W3-480.txt')
FFA_36 = np.loadtxt('FFA-W3-360.txt')
FFA_301 = np.loadtxt('FFA-W3-301.txt')
FFA_241 = np.loadtxt('FFA-W3-241.txt')
static_blade_deflection_data = np.loadtxt('bladestruc.txt')
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
EI_1 = np.zeros(len(static_blade_deflection_data))
EI_2 = np.zeros(len(static_blade_deflection_data))
nu = np.zeros(len(static_blade_deflection_data)) #structural pitch
mass = np.zeros(len(static_blade_deflection_data))

for i in range(len(blade_data)):
    r_data[i] = blade_data[i][0]
    c_data[i] = blade_data[i][1]
    b_data[i] = blade_data[i][2]
    tc_data[i] = blade_data[i][3]

for i in range(len(static_blade_deflection_data)):
    EI_1[i] =  static_blade_deflection_data[i][3]
    EI_2[i] =  static_blade_deflection_data[i][4]
    nu[i] =  static_blade_deflection_data[i][1]
    mass[i] = static_blade_deflection_data[i][2]


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

r_new = np.linspace(blade_data[0][0],R,40)
c_new = np.interp(r_new,r_data,c_data)
b_new = np.interp(r_new,r_data,b_data)
tc_new = np.interp(r_new,r_data,tc_data)
EI_1_new = np.interp(r_new,r_data,EI_1)
EI_2_new = np.interp(r_new,r_data,EI_2)
nu_new = np.interp(r_new,r_data,nu)
mass_new = np.interp(r_new,r_data,mass)

B = 3
rho = 1.225
V_0 = 6.0
dx = 1e-5
Theta_p = -0.68*np.pi/180 #deg to rad
#Theta_p = 0*0.896*np.pi/180
#-----------------------BEM algorithm-----------------------------------------    
omega = 0.5134243994347621
#omega=0.9253
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
                
                aoa = phi - (Theta_p + beta)
                
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
    
                if abs(a_old - a) < 1e-6 and abs(a_prime_old - a_prime) < 1e-6:
                    
                    k = 1
            #print(C_l)
            aoa = phi - (Theta_p + beta)
        
            V_rel_sq = (V_0*(1-a))**2 + (omega*r*(1+a_prime))**2
            
            C_n = C_l*np.cos(phi) + C_d*np.sin(phi)

            C_t = C_l*np.sin(phi) - C_d*np.cos(phi)

            p_n[i] = (1/2)*rho*V_rel_sq*chord*C_n

            p_t[i] = (1/2)*rho*V_rel_sq*chord*C_t
            
            F = (2/np.pi)*np.arccos(np.exp(-(B*(R-r))/(2*r*np.sin(abs(phi)))))
        
integral = trapezoid(r_new*p_t, r_new, dx)
integral_T = trapezoid(p_n, r_new, dx)
P = omega*B*integral
Thrust = B*integral_T
#print('p_t:',p_t)
#print('p_n:',p_n)
#-----------------------------------Ashes data--------------------------------------------------------------------------
p_n_Ashes = np.array([56.083, 67.4511, 161.475, 185.261, 203.344, 231.522, 275.852, 522.255, 1189.12, 1351.28, 1339.67, 1473.41, 1605.27, 1740.89, 1878.84, 2017.26, 2158.7, 2298.05, 2429.98, 2558.77, 2678.17, 2785.24, 2878.42, 2956.64, 3025.93, 3079.95, 3129.22, 3158.94, 3185.42, 3192.96, 3189.44, 3167.22, 3123.58, 3048.09, 2938.38, 2800.05, 2596.37, 2302.41, 1858.65, 0])
p_t_Ashes = np.array([-13.421, -31.3817, 70.7997, 59.27, 42.5693, 30.495, 22.6857, 77.776, 290.413, 289.744, 296.335, 296.035, 295.355, 298.892, 298.34, 297.691, 296.949, 296.097, 295.167, 294.147, 293.07, 291.937, 290.696, 289.239, 287.378, 284.988, 281.888, 277.919, 272.918, 266.63, 258.929, 249.594, 238.371, 224.926, 208.875, 189.621, 166.148, 136.504, 95.5187, 0])
r_Ashes = np.array([0, 2.64302, 5.37977, 8.20274, 11.1031, 14.071, 17.0953, 20.1641, 23.2647, 26.3837, 29.5076, 32.6228, 35.7156, 38.773, 41.7824, 44.732, 47.6111, 50.4099, 53.1201, 55.7344, 58.247, 60.6534, 62.9501, 65.1352, 67.2076, 69.1675, 71.0159, 72.7545, 74.386, 75.9133, 77.3402, 78.6705, 79.9085, 81.0585, 82.1252, 83.113, 84.0265, 84.8703, 85.6487, 86.366])
#------------------------------------plotting loads---------------------------------------------------------------------------
fig, ax1 = plt.subplots()
ax1.plot(r_Ashes + r_new[0],p_t_Ashes,'r-o',label='Ashes')
ax1.plot(r_new,p_t, 'b-o', label='BEM')
ax1.set_title(r"$p_t$ distribution comparison for $V_0 = 10.0$ m/s")
ax1.tick_params(axis='both', labelsize=16)  
ax1.set_xlabel(r"$ r $ [m]", fontsize=16)
ax1.set_ylabel(r"$ p_t $[N/m]", fontsize=16)
ax1.legend()
ax1.grid(True, alpha=0.3)

fig, ax2 = plt.subplots()
ax2.plot(r_Ashes + r_new[0], p_n_Ashes,'r-o',label='Ashes')
ax2.set_title(r"$p_n$ distribution comparison for $V_0 = 10.0$ m/s")
ax2.plot(r_new,p_n,'b-o',label='BEM')
ax2.tick_params(axis='both', labelsize=16)  
ax2.set_xlabel(r"$ r $ [m]", fontsize=16)
ax2.set_ylabel(r"$ p_n $[N/m]", fontsize=16)
ax2.legend()
ax2.grid(True, alpha=0.3)
#---------------------------------------------calculating deflections-------------------------------------------------------
da_y = deflections(p_n, p_t, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[1]
da_z = deflections(p_n, p_t, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[2]
d_y = deflections(p_n, p_t, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[3]
d_z = deflections(p_n, p_t, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[4]
data_out = np.column_stack((r_new, da_y, da_z, d_y, d_z))
np.savetxt("deflections_output.txt", data_out, header="r[m]   da_y   da_z   d_y   d_z", fmt="%.6f")
data_out_loads = np.column_stack((r_new, p_n, p_t))
np.savetxt("loads.txt", data_out_loads, header="r[m]   p_n   p_t", fmt="%.6f")
#print(D)
plt.figure()
plt.plot(r_new, da_y, '-o')
plt.xlabel('r [m]')
plt.ylabel('da_y')
plt.title('deflection angle da_y along blade (root -> tip)')
plt.grid(True)

plt.figure()
plt.plot(r_new, da_z, '-o')
plt.xlabel('r [m]')
plt.ylabel('da_z')
plt.title('deflection angle da_y along blade (root -> tip)')
plt.grid(True)

plt.figure()
plt.plot(r_new, d_y, 'r-o', label='d_y')
plt.plot(r_new, d_z, 'b-o', label='d_z')
plt.xlabel('r [m]')
plt.ylabel('d')
plt.title('deflection d_y along blade (root -> tip)')
plt.legend()
plt.grid(True)

#plt.figure()
#plt.plot(r_new, d_z, '-o')
#plt.xlabel('r [m]')
#plt.ylabel('d_z')
#plt.title('deflection d_z along blade (root -> tip)')
#plt.grid(True)
#---------------------------------------------eigenmode calculation---------------------------------------------------------
ni = len(r_new)             # total grid points (same as MATLAB)
nd = ni - 1                 # number of active DOFs per direction
DOF = 2 * nd                # total DOFs

# Initialize F
F = np.zeros((DOF, DOF))

# ========== FIRST LOOP: load in y-direction ==========
for j in range(1, ni):                 # MATLAB: j = 2:ni
    p_y = np.zeros(ni)
    p_z = np.zeros(ni)
    p_y[j] = 1.0                       # load at node j (Python index matches)

    d = deflections(p_z, p_y, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[0]
    d = np.ravel(d)                    # length = 2*(ni-1)
    
    # Fill column j-1
    F[0:nd, j-1]   = d[0:nd]           # uy(2:end)
    F[nd:, j-1]    = d[nd:]            # uz(2:end)
    #if F[j, ] @ p_y == d:
     #   print("everything fine")
# ========== SECOND LOOP: load in z-direction ==========
for j in range(1, ni):
    p_y = np.zeros(ni)
    p_z = np.zeros(ni)
    p_z[j] = 1.0

    d = deflections(p_z, p_y, r_new, b_new, nu_new, Theta_p, EI_1_new, EI_2_new)[0]
    d = np.ravel(d)

    col = nd + j - 1                   # MATLAB: ni + j - 2

    F[0:nd, col] = d[0:nd]
    F[nd:, col]  = d[nd:]
    
M = np.zeros((DOF, DOF))

for i in range(1, ni):   # MATLAB: i = 2:ni
    M[i-1, i-1] = mass_new[i]
    M[nd + i - 1, nd + i - 1] = mass_new[i]

condF = np.linalg.cond(F)
print("cond(F) =", condF)
A = F@M
print("Any NaN in A?:", np.isnan(A).any())
print("Any Inf in A?:", np.isinf(A).any())
print("Rank of F:", np.linalg.matrix_rank(F))
print("Size of F:", F.shape[0])
print("NaN in F:", np.isnan(F).any())
print("NaN in M:", np.isnan(M).any())
print("NaN in F @ M:", np.isnan(F @ M).any())
#print(M)
eigenvalue,eigenvector = eig(A)
#print('eigenvalue',eigenvalue)
#print('eigenvector',eigenvector[:,3])
eigenfrequency = np.sqrt(1 / eigenvalue)
print('eigenfrequency',eigenfrequency)
mode_shape_1y = eigenvector[:, 0][0:nd]/np.max(np.abs(eigenvector[:, 0][0:DOF]))
mode_shape_1z = eigenvector[:, 0][nd:]/np.max(np.abs(eigenvector[:, 0][0:DOF]))
plt.figure()
plt.plot(r_new[1:], mode_shape_1y, 'r-o', label='d_y')
plt.plot(r_new[1:], mode_shape_1z, 'b-o', label='d_z')
plt.xlabel('r [m]')
plt.ylabel('d/dmax')
plt.title('deflection d along blade (mode 1) (root -> tip)')
plt.legend()
plt.grid(True)
mode_shape_2y = eigenvector[:, 1][0:nd]/np.max(np.abs(eigenvector[:, 1][0:DOF]))
mode_shape_2z = eigenvector[:, 1][nd:]/np.max(np.abs(eigenvector[:, 1][0:DOF]))
plt.figure()
plt.plot(r_new[1:], mode_shape_2y, 'r-o', label='d_y')
plt.plot(r_new[1:], mode_shape_2z, 'b-o', label='d_z')
plt.xlabel('r [m]')
plt.ylabel('d/dmax')
plt.title('deflection d along blade (mode 2) (root -> tip)')
plt.legend()
plt.grid(True)
mode_shape_3y = eigenvector[:, 2][0:nd]/np.max(np.abs(eigenvector[:, 2][0:DOF]))
mode_shape_3z = eigenvector[:, 2][nd:]/np.max(np.abs(eigenvector[:, 2][0:DOF]))
plt.figure()
plt.plot(r_new[1:], mode_shape_3y, 'r-o', label='d_y')
plt.plot(r_new[1:], mode_shape_3z, 'b-o', label='d_z')
plt.xlabel('r [m]')
plt.ylabel('d/dmax')
plt.title('deflection d along blade (mode 3) (root -> tip)')
plt.legend()
plt.grid(True)

