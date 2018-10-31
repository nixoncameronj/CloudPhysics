
# coding: utf-8

# In[11]:


import numpy as np
import Bolton

import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from readsoundings import parse_SPC

sounding = '/home/cjnix/cjnix/Documents/Python/git_repos/CloudPhysics/april9sounding'
sounding_data = parse_SPC(sounding)



# In[18]:


C_to_K = 273.15

# serves as alpha
skew_slope = 40

def x_from_Tp(T, p):
    '''calculates x(T,p)
    '''
    x = T - (skew_slope * np.log(p))
    return x

def y_from_p(p):
    '''calculates y(p)
    '''
    y = -(np.log(p))
    return y

def p_from_y(y):
    '''calculates p(y)
    '''
    p = np.exp(-y)
    return p

def T_from_xp(x, p):
    '''calculates T(x,p)
    '''
    T = x + (skew_slope * np.log(p))
    return T

def to_thermo(x, y):
    '''Transform (x,y) coordinates to T in degrees Celsius
       and p in mb'''
    p = p_from_y(y)
    T_C = T_from_xp(x, p) - C_to_K
    return T_C, p

def from_thermo(T_C, p):
    '''Transform T_C (in degrees Celsius)
       and p (in mb) to (x,y)'''
    y = y_from_p(p)
    x = x_from_Tp(T_C+C_to_K, p)
    return x, y

# def theta_e(T, p, p_0):
#     '''Calculates theta_e from T, p, p_0
#     '''
#     R_d = 287.1
#     c_l = 4218
#     c_p = 1005
#     w_s = Bolton.sat_mixing_ratio(p, T)
#     w_t = w_s
#     C_wd = c_pd + (w_s*c_l)
#     L_v = (3.139*(10**6)) - (c_l - c_p)*T
#     theta_e = T*((p/p_0)**(R_d/c_wd))*np.exp((L_v*w_s)/(c_wd*T))
#     return theta_e

# def theta_e_lines(T, p, p_0=1000.0):
#     '''calculates theta-e lines
#     '''
#     w = sat_mixing_ratio(p, T)
#     theta_e_lines = theta_e(T, p, w, p_0)
#     return theta_e_lines

# values along the bottom and left edges
p_bottom = 1050.0
p_top = 150
T_min = -40 + C_to_K
T_max = 50 + C_to_K

# calculate graph bounds
x_min = x_from_Tp(T_min,p_bottom)
x_max = x_from_Tp(T_max,p_bottom)
y_min = y_from_p(p_bottom)
y_max = y_from_p(p_top)

# define constant levels
p_levels = np.arange(1000, 150-50, -50)
T_C_levels = np.arange(-120, 40+10, 10)
T_levels = T_C_levels + C_to_K
P_levels = np.arange(1000, 150-50, -50)
theta_levels = (np.arange(-40, 100+10, 10)) + C_to_K
theta_e_levels = (np.arange(-40, 100+10, 10)) + C_to_K
theta_ep_levels = theta_levels.copy()
mixing_ratios = np.asarray([.0004,.001,.002,.003,.005,.008,.012,.016,.020])


# In[19]:


p_all = np.arange(p_bottom,p_top-1,-1)
y_p_levels = y_from_p(p_levels)
y_all_p = y_from_p(p_all)
x_T_levels = [x_from_Tp(Ti, p_all) for Ti in T_levels]
x_thetas = [x_from_Tp(Bolton.theta_dry(theta_i, p_all),p_all) for theta_i in theta_levels]
x_mixing_ratios = [x_from_Tp(Bolton.mixing_ratio_line(p_all,mixing_ratio_i)+C_to_K,p_all) for mixing_ratio_i in mixing_ratios]
x_thetaes = [x_from_Tp(Bolton.theta_dry(thetae_i, p_all),p_all) for thetae_i in theta_e_levels]

mesh_T, mesh_p = np.meshgrid(np.arange(-60.0, T_levels.max()-C_to_K+0.1, 0.1), p_all)
theta_ep_mesh = Bolton.theta_ep_field(mesh_T, mesh_p)
theta_e_mesh = Bolton.theta_e_field(mesh_T, mesh_p)

# sounding data
snd_p = sounding_data['p']
good_p = (snd_p > 200) & (snd_p < 1000)
y_snd_p = y_from_p(snd_p)

snd_T = sounding_data['T']
# all temperature values, deg. C, should be in this range.
good_T = (snd_T > -100.0) & (snd_T < 60.0)
x_snd_T = x_from_Tp(snd_T, snd_p)+C_to_K
y_snd_T = y_from_p(snd_p)

snd_Td = sounding_data['Td']
good_Td = (snd_Td > -100.0) & (snd_Td < 60.0)
x_snd_Td = x_from_Tp(snd_Td, snd_p)+C_to_K
y_snd_Td = y_from_p(snd_p)

#


# In[20]:


%matplotlib inline
skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))
fig = plt.figure()
ax = Subplot(fig,1,1,1,grid_helper = skew_grid_helper)
def format_coord(x, y):
    T, p = to_thermo(x, y)
    return "{0:5.1f} C, {1:5.1f} mb".format(float(T),float(p))
ax.format_coord = format_coord
fig.add_subplot(ax)

for yi in y_p_levels:
    ax.plot((x_min, x_max), (yi,yi), color=(1.0, 0.8, 0.8))

for x_T in x_T_levels:
    ax.plot(x_T, y_all_p, color=(1.0, 0.5, 0.5))

#doesn't work!
for x_theta in x_thetas:
    ax.plot(x_theta, y_all_p, color=(1.0, 0.7, 0.7))

for x_mixing_ratio in x_mixing_ratios:
    good = p_all >= 600 # restrict mixing ratio lines to below 600 mb
    ax.plot(x_mixing_ratio[good], y_all_p[good], color=(0.8, .8, 0.6))

n_moist = len(theta_ep_levels)
moist_colors = ((0.6,0.9,0.7),)*n_moist
ax.contour(x_from_Tp(mesh_T+C_to_K, mesh_p), y_from_p(mesh_p),
    theta_ep_mesh, theta_ep_levels, colors=moist_colors)

#doesn't work!
theta_e_lines = len(theta_e_levels)
moist_colors = ((0.6,0.9,0.7),)*theta_e_lines
ax.contour(x_from_Tp(mesh_T+C_to_K, mesh_p), y_from_p(mesh_p),
    theta_e_mesh, theta_e_levels, colors='lightblue')

ax.plot(x_snd_Td, y_snd_p, linewidth=2, color='g')
ax.plot(x_snd_T, y_snd_p, linewidth=2, color='r')
# your code for plotting theta_e (reversible)

ax.axis((x_min, x_max, y_min, y_max))

#ax.set_xlim(-40,0)

plt.show()
