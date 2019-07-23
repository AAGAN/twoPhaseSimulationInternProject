import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class container:
    """container class to calculate the properties of a pressurized container"""

    def __init__(self,
                agentInitialTemp, #initial temperature of the gas inside the cylinder
                initialMass,
                cylInitialTemp,
                cylLength,
                cylDiam,
                nozzleDiam,
                wallThickness,
                ambientTemp,
                backPressure,
                endTime,
                numTimeSteps = 1000):
        #define the parameters
        self.T_a0 = agentInitialTemp 
        self.m_a0 = initialMass #initial mass of gas inside the cylinder
        self.T_w0 = cylInitialTemp #initial temperature of the cylinder walls
        self.L = cylLength #length of the cylinder
        self.D = cylDiam #diameter of the cylinder
        self.A_s = np.pi * self.D * ( self.L + self.D / 2.0 ) #surface area of the walls of the tank
        self.D_e = nozzleDiam #diameter of the nozzle at the exit (assuming a converging nozzle here)
        self.A_e = np.pi * self.D_e **2 / 4.0 #exit area of the nozzle
        self.A_c = np.pi* self.D ** 2 /4.0 #Area of the cylinder
        self.V = self.A_c * self.L # volume inside the tank 
        self.R_a = 287 #gas constant for gas
        self.c_w = 1250 #specific heat of the tank wall material
        self.g_c = 1 #proportionality appearing in Newton's second law 
        self.g = 9.81 #gravitational acceleration
        self.J = 1 #joule's constant
        self.rho_container = 1200 #container wall's density 
        self.wall_thickness = wallThickness #container wall's thickness
        self.m_w = self.rho_container * (np.pi * (self.D+self.wall_thickness/2.0)*self.L*self.wall_thickness + 2 * np.pi / 4.0 * self.D ** 2 * self.wall_thickness) #mass of cylinder
        self.T_inf = ambientTemp #ambient temperature
        self.p_b = backPressure #back pressure (ambient pressure)
        self._endTime = endTime
        #self.tspan = np.linspace(0,endTime,numTimeSteps) #

    def mu(self, T):#=============================================function of film temperature
        return 1.8e-3 #* u.kilogram / (u.meter * u.second)

    def k_t(self, T):#=============================================function of film temperature
        return 0.026 #* u.watt / u.meter / u.kelvin #thermal conductivity of the gas==============================

    def beta(self, T):
        return 1/T #for ideal gasses 

    def p_a(self, m_a,R_a,T_a,V):
        p = m_a * R_a * T_a / V
        return p

    def c_p(self, T):#=============================================function of film temperature
        return 1000 #* u.joule / u.kilogram / u.kelvin

    def c_v(self, T):#=============================================function of film temperature
        return 718 #* u.joule / u.kilogram / u.kelvin

    def k(self, T):
        return self.c_p(T)/self.c_v(T)
        
    def hbar_i(self, T_w, T, ma):
        _beta = self.beta(T)
        _c_p = self.c_p((T+T_w)/2.0)
        _rho = ma / self.V
        _mu = self.mu((T+T_w)/2.0)
        _k_t = self.k_t((T+T_w)/2.0)
        _neu = _mu / _rho
        Gr = (self.g * _beta * abs(T_w - T) * self.L ** 3 ) / _neu ** 2
        Pr = _c_p * _mu / _k_t
        #print(Gr,'\n',_mu,'\n',Pr,'\n',(Gr*Pr))
        if 10**4 < Gr * Pr <= 10**9:
            C = 0.59
            m = 0.25
        else:
            C = 0.129
            m = 1.0/3.0
        
        h = (C * (Gr * Pr) ** m)*_k_t/self.L
        return h

    def hbar_o(self, T_w):
        _beta = self.beta(self.T_inf)
        _c_p = self.c_p((self.T_inf+T_w)/2.0)
        _rho = self.p_b / self.R_a / self.T_inf
        _mu = self.mu((self.T_inf+T_w)/2.0)
        _neu = _mu / _rho
        _k_t_inf = self.k_t((self.T_inf+T_w)/2.0)
        Gr = (self.g * _beta * abs(T_w - self.T_inf) * self.L ** 3 ) / _neu ** 2
        Pr = _c_p * _mu / _k_t_inf
        #print(Gr,'\n',_mu,'\n',Pr,'\n',(Gr*Pr))
        if 10**4 < Gr * Pr <= 10**9:
            C = 0.59
            m = 0.25
        else:
            C = 0.129
            m = 1.0/3.0
        
        h = (C * (Gr * Pr) ** m)*_k_t_inf/self.L
        return h

    def mdot_a(self, m_a, T_a, T_w):
        _k = self.k(T_a)
        _p_a = self.p_a(m_a,self.R_a,T_a,self.V)
        if self.p_b/_p_a > 0.528:
            p_e = self.p_b
            Me = 2 / (_k-1) * (1-(p_e/_p_a)**((_k-1)/_k))
            Te = T_a / (1+(_k-1)/2*Me**2)
            ce = (_k*self.R_a*self.g_c*Te)**0.5
            ve = Me * ce
            mdote = -p_e/self.R_a/Te*self.A_e*ve
        else: #chocked flor for p_b/p_a <= 0.528
            p_e = 0.528 * _p_a
            Me = 1.0
            Te = T_a / (1+(_k-1)/2.0)
            ve = ce = (_k*self.R_a*self.g_c*Te)**0.5
            mdote = -p_e/self.R_a/Te*self.A_e*ve
        return mdote

    def T_e(self, m_a, T_a, T_w):
        _k = self.k(T_a)
        _p_a = self.p_a(m_a,self.R_a,T_a,self.V)
        if self.p_b/_p_a > 0.528:
            p_e = self.p_b
            Me = 2 / (_k-1) * (1-(p_e/_p_a)**((_k-1)/_k))
            Te = T_a / (1+(_k-1)/2*Me**2)
        else: #chocked flor for p_b/p_a <= 0.528
            p_e = 0.528 * _p_a
            Me = 1.0
            Te = T_a / (1+(_k-1)/2.0)
        return Te

    def v_e(self, m_a, T_a, T_w):
        _k = self.k(T_a)
        _p_a = self.p_a(m_a,self.R_a,T_a,self.V)
        if self.p_b/_p_a > 0.528:
            p_e = self.p_b
            Me = 2 / (_k-1) * (1-(p_e/_p_a)**((_k-1)/_k))
            Te = T_a / (1+(_k-1)/2*Me**2)
            ce = (_k*self.R_a*self.g_c*Te)**0.5
            ve = Me * ce
        else: #chocked flor for p_b/p_a <= 0.528
            p_e = 0.528 * _p_a
            Me = 1.0
            Te = T_a / (1+(_k-1)/2.0)
            ve = ce = (_k*self.R_a*self.g_c*Te)**0.5
        return ve

    def dTa_dt(self, m_a, T_a, T_w):
        _hbar_i = self.hbar_i(T_w, T_a, m_a)
        _mdot_e = -self.mdot_a(m_a, T_a, T_w)
        _T_e = self.T_e(m_a, T_a, T_w)
        _v_e = self.v_e(m_a, T_a, T_w)
        firstTerm = _mdot_e / m_a * T_a
        secondTerm = self.A_s * _hbar_i / self.c_v(T_a) / m_a * (T_w - T_a)
        thirdTerm = -_mdot_e / m_a * self.c_p(T_a) / self.c_v(T_a) * _T_e
        forthTerm = -_mdot_e * _v_e ** 2 / 2.0 / self.J / self.g_c / self.c_v(T_a) / m_a
        dTa_dt = firstTerm + secondTerm + thirdTerm + forthTerm
        return dTa_dt # firstTerm , secondTerm , thirdTerm , forthTerm

    def dTw_dt(self, m_a,T_a,T_w):
        _hbar_o = self.hbar_o(T_w)
        _hbar_i = self.hbar_i(T_w, T_a, m_a)
        firstTerm = self.A_s*_hbar_o/self.m_w/self.c_w*(self.T_inf-T_w) 
        secondTerm = -self.A_s*_hbar_i/self.m_w/self.c_w*(T_w-T_a)
        dTw_dt = firstTerm+secondTerm
        #return dTw_dt
        return dTw_dt

    #RK4 implementation
    def f(self,t, y):
        dydt = [self.mdot_a(y[0], y[1], y[2]), self.dTa_dt(y[0], y[1], y[2]) , self.dTw_dt(y[0], y[1], y[2])]
        return dydt

    def solve(self):
        yinit = [self.m_a0,self.T_a0,self.T_w0]
        self.sol = solve_ivp(lambda t, y: self.f(t, y), 
                         [0, self._endTime], yinit,method = 'Radau', jac = None)
        # self.sol = solve_ivp(lambda t, y: self.f(t, y), 
        #                 [self.tspan[0], self.tspan[-1]], yinit,method = 'Radau', jac = None)

    def t(self):
        return self.sol.t

    def ma(self):
        return self.sol.y[0]

    def Ta(self):
        return self.sol.y[1]

    def Tw(self):
        return self.sol.y[2]

    def Pt(self):
        return self.p_a(self.sol.y[0],self.R_a,self.sol.y[1],self.V)

    def Ve(self):
        ve = np.zeros(len(self.sol.y[0]))
        for i in range(len(self.sol.y[0])):
            ve[i] = self.v_e(self.sol.y[0][i],self.sol.y[1][i],self.sol.y[2][i])
        return ve

    def Me(self):
        Me = np.zeros(len(self.sol.y[0]))
        for i in range(len(self.sol.y[0])):
            Me[i] = self.M_e(self.sol.y[0][i],self.sol.y[1][i],self.sol.y[2][i])
        return Me

    def plot(self):
        plt.rcParams['figure.figsize'] = (15, 10)
        plt.rcParams.update({'font.size': 15})
        plt.subplot(2,2,1)
        plt.plot(self.sol.t, self.sol.y[0],linewidth = 3)
        plt.ylabel("$m_a (kg)$")
        plt.ylim(0)
        plt.grid()
        plt.subplot(2,2,2)
        plt.ylabel("$T_a (K)$")
        plt.plot(self.sol.t, self.sol.y[1],linewidth = 3)
        plt.grid()
        #plt.legend()
        #plt.xlabel("time (s)")
        plt.subplot(2,2,3)
        plt.ylabel("$T_w (K)$")
        plt.plot(self.sol.t, self.sol.y[2],linewidth = 3)
        plt.grid()
        #plt.legend()
        plt.xlabel("time (s)")
        plt.subplot(2,2,4)
        plt.ylabel("$Pressure (psi)$")
        plt.plot(self.sol.t, self.p_a(self.sol.y[0],self.R_a,self.sol.y[1],self.V)*0.000145038,linewidth = 3)
        plt.grid()
        #plt.legend()
        plt.xlabel("time (s)")
        plt.show()

    def M_e(self, m_a, T_a, T_w):
        _k = self.k(T_a)
        _p_a = self.p_a(m_a,self.R_a,T_a,self.V)
        if self.p_b/_p_a > 0.528:
            p_e = self.p_b
            Me = 2 / (_k-1) * (1-(p_e/_p_a)**((_k-1)/_k))
        else: #chocked flor for p_b/p_a <= 0.528
            p_e = 0.528 * _p_a
            Me = 1.0
        return Me

    def plotVe(self):    
        ve = np.zeros(len(self.sol.y[0]))
        for i in range(len(self.sol.y[0])):
            ve[i] = self.M_e(self.sol.y[0][i],self.sol.y[1][i],self.sol.y[2][i])
        plt.plot(self.sol.t, ve,linewidth = 3)
        plt.grid()
        #plt.legend()
        plt.show()
