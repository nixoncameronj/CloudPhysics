
# coding: utf-8

# In[4]:


import numpy as np


# In[7]:


C_to_K = 273.15
c_p_dry = 1005.7
c_V_dry = 4190
eps = 0.6220
k_dry = 0.2854

def sat_vapor_pressure(T):
    '''calculates saturation vapor pressure from T
    '''
    es = 6.112*np.exp((17.67*T)/(T+243.5))
    return es

def sat_vapor_temperature(e_s):
    '''calculates saturation vapor temperature from e_s
    '''
    T = (243.5*np.log(e_s) - 440.8)/(19.48-np.log(e_s))
    return T

def sat_mixing_ratio(p,T):
    '''calculates saturation mixing ratio from p, T
    '''
    e_s = sat_vapor_pressure(T)
    w = eps*(e_s/(p-e_s))
    return w

def mixing_ratio_line(p, w_s):
    '''calculates temperature from p, w_s
    '''
    e_s = p*((w_s/eps)+1)
    T = sat_vapor_temperature(e_s)
    return T

def RH(T, p, w):
    '''calculates relative humidity from T, p, w
    '''
    w_s = sat_mixing_ratio(p,T)
    RH = (w/w_s)*100
    return RH

def T_LCL(T, RH):
    '''calculates LCL temperature in Kelvin
    '''
    T_LCL = (1/((1/(T+C_to_K-55))-(np.log(RH/100)/2840)))+55
    return T_LCL

# accept many pressures?
def theta_dry(theta, p, p_0=1000.0):
    '''calculates theta dry from theta, p and reference pressure
    '''
    theta_dry = theta*((p_0/(p-0))**k_dry)
    return theta_dry + C_to_K

def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    '''calculates theta ep from T, p, w and reference pressure
    '''
    rh = RH(T, p, w)
    t_lcl = T_LCL(T, rh)
    pseudoeq_potential_T = ((T+C_to_K)*((p_0/p)**(0.2854*(1-0.28*np.exp(-3*w)))))*(np.exp((3.376/t_lcl)-0.00254))*(w*(1+(w*0.81*np.exp(-3))))
    return pseudoeq_potential_T

# probably not right
def theta_ep_field(T, p, p_0=1000.0):
    '''calculates moist adiabats
    '''
    w = sat_mixing_ratio(p, T)
    theta_ep_field = pseudoeq_potential_T(T, p, w) + C_to_K
    return theta_ep_field

