import numpy as np
from scipy.optimize import fsolve
## Function for Orifice 

#For forward pass (Things already known prior to calculation is the mdot, M1, P01, T01 i.e upstream parameters)
def mass_critical (P1,T1,A_o,y,Cd,Mw_gas):
    P_1=P1
    T_1=T1
    _A_o=A_o
    C_d=Cd
    _y=y
    R= 8.3145
    _Mw_gas=Mw_gas
    _R=R/Mw_gas
    m_c= C_d*_A_o*P_1*(2/_R/T_1)**0.5*(_y/(_y+1)*(2/(_y+1))**(2/(_y-1)))**0.5
    return m_c

def stat_P (P0, M, y):
    P_0=P0
    P_stat= P_0/(1+(y-1)/2*M**2)**(y/(y-1))
    return P_stat

def stat_T(T0, M, y):
    T_0=T0
    T_stat= T_0/(1+(y-1)/2 *M**2)
    return T_stat

def main_nozzle_forward(P01, T01, mdot,M, d_orifice, d_upstream_pipe,y, Mw_gas):
    A_o= (np.pi*d_orifice**2/4) #Area of orifice
    A_p= (np.pi*d_upstream_pipe**2/4) #Area of pipe 

    Cd= 0.61 #Coefficient of discharge, It will be less than 1
    #y= 1.4  gamma value changes with the type of gas and conditions
    #g=9.81 #gravitational constant
    R= 8.3145
    _Mw_gas = Mw_gas #since input in grams, if in kg remove this
    _R=R/_Mw_gas
    B = np.sqrt(A_o/A_p) 
    
    P_01=P01
    T_01=T01
    _m=mdot
    
    P_1=stat_P(P_01,M,y) #these equatiions need Mach number
    T_1=stat_T(T_01,M,y)
    
    rho1= P_1/(_R*T_1)
    
    m_c= mass_critical(P_1, T_1, A_o,y,Cd,Mw_gas)
    
    if _m>m_c:
        r_f = m_c/_m #reduction factor
        print("mass flow rate is greater than critical mass flow rate")
        print ('the output you got is the ratio needed to be multiplied to the original first m assumed') 
        return r_f
    elif _m == m_c:
        print ('Flow is choked')
        P_2 = P_1*(2/(y+1))**(y/(y-1))
        P_02 = P_01 + P_2 - P_1 # need to calculate this? or P_02 is P_2 critical itself?
        return P_2,P_02
    else: #no choking
        A=A_o
        _B=B
        def equation(P_2):
            return (Cd*(1-(0.333+1.145*(_B**2+0.7*_B**5+ 12*_B**13))*(P_1-P_2)/(y*P_1))*A*(2*rho1*(P_1-P_2))**0.5-_m)

        _P_2= fsolve(equation, P_1-2000 , xtol=0.000001)
        #print("Y",(1-(0.333+1.145*(_B**2+0.7*_B**5+ 12*_B**13))*(P_1-_P_2)/(y*P_1)))
        #print ('P2=',_P_2,'error=', equation(_P_2))
        P_02 = P_01 + _P_2 - P_1
        return _P_2[0],P_02[0] # we can even get M and P2 from here
    

    
