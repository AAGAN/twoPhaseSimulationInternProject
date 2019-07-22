import numpy as np
from scipy.integrate import solve_ivp

class calcTw:
    def __init__(self,
                L, #length of the pipe
                D, #internal Diameter of the pipe
                M1, #Initial Mach number
                ff, #friction factor of the pipe
                T01, #initial total temperature of the gas
                T1, #initial temperature of the gas
                P01, #initial pressure of the gas
                P1, #initial pressure of the gas
                rho1, #initial density of the gas
                V1, #initial velocity of the gas
                gamma, #ratio of specific heats
                Tw, #wall temperature for the case of constant wall temperature
                numSteps = 10000 #number of steps from the begining to the end of the pipe
    ):
        self._L = L
        self._D = D
        self._M1 = M1
        self._ff = ff
        self._T01 = T01
        self._T1 = T1
        self._P01 = P01
        self._P1 = P1
        self._gamma = gamma
        self._Tw = Tw
        self._numSteps = numSteps
        self._rho1 = rho1
        self._V1 = V1

    def FuncT0(self, Msq, Gamma):
        FuncT0 = Msq*(1.0+Gamma*Msq)*(1.0+(Gamma-1.0)/2.0*Msq)/(1.0-Msq)
        return FuncT0

    def FuncFf(self, Msq, Gamma):
        FuncFf = Gamma*Msq**2*(1.0+(Gamma-1.0)/2.0*Msq)/(1.0-Msq)
        return FuncFf

    def dT0_dx(self, T0):
        _dT0_dx = 2 * self._ff / self._D * (self._Tw - T0)
        return _dT0_dx

    def dM2_dx(self, Msq, T0):
        _dM2_dx = self.FuncT0(Msq, self._gamma) * self.dT0_dx(T0) / T0 + self.FuncFf(Msq, self._gamma) * 4.0 * self._ff / self._D
        return _dM2_dx

    def f(self, t,y):
        dydt = [self.dM2_dx(y[0], y[1]),self.dT0_dx(y[1])]
        return dydt

    def solve(self):
        yinit = [self._M1 **2, self._T01]
        # tspan = np.linspace(0,self._L,self._numSteps)
        # self._soln = solve_ivp(lambda t, y: self.f(t, y), [tspan[0], tspan[-1]], yinit, method = 'Radau', t_eval=tspan, jac = None)
        self._soln = solve_ivp(lambda t, y: self.f(t, y), [0, self._L], yinit, method = 'Radau', jac = None)
        self._M2 = np.sqrt(self._soln.y[0][-1])
        self._T02 = self._soln.y[1][-1]
        self.T2()
        self.P2()
        self.P02()
        self.rho2()
        self.V2()

    def T2(self):
        self._T2 = self._T02 / (1+(self._gamma-1.0)/2.0*self._M2**2)

    def P2(self):
        self._P2 = self._P1 * self._M1/self._M2*np.sqrt((1+(self._gamma-1)/2*self._M1**2)/(1+(self._gamma-1)/2*self._M2**2))*np.sqrt(self._T02/self._T01)

    def P02(self):
        comg = (self._gamma - 1 ) / 2
        self._P02 = self._P01 * self._P2 / self._P1 * ((1+comg*self._M2**2)/(1+comg*self._M1**2))**(self._gamma/(self._gamma-1))

    def rho2(self):
        self._rho2 = self._rho1 * self._P2/self._P1 * self._T1 / self._T2

    def V2(self):
        self._V2 = self._M2 / self._M1 * np.sqrt(self._T2/self._T1)