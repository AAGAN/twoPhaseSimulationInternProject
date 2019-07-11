import pint
import numpy as np
import matplotlib.pyplot as plt
u = pint.UnitRegistry()
from scipy.integrate import solve_ivp

class pipe:
    def __init__(self,
                pipeLength,
                initialMach,
                wallTemp,
                gamma):
        self.L = pipeLength
        self.M1 = initialMach
        self.gamma = gamma
        self.h = 0.01
        self.Tw01 = wallTemp
        self.Eps = 1e-6 #convergence criteria

    def findT02(self):
        OrM1 = self.M1
        M2 = self.M1 + self.h
        Prod = 1.0
        CommG = (self.gamma - 1.0)/2.0
        while (M2<1):
            Soln = self.Solve(self.M1,M2,self.Tw01,self.gamma,self.Eps)
            if Soln>1.0:
                Prod = Prod * Soln
                X = 2.0 * (np.log10(self.Tw01-1.0)-np.log10(self.Tw01-Prod))
                T = Prod * (1.0+CommG*OrM1*OrM1)/(1.0+CommG*M2*M2)
                P = OrM1/M2*T**0.5
                P0 = Prod**0.5*OrM1/M2*((1.0+CommG*M2*M2)/(1.0+CommG*OrM1*OrM1))**((self.gamma+1.0)/(2.0*(self.gamma-1.0)))
                print(M2,Prod,X,T,P,P0)
            self.M1 = M2
            M2 = M2 + self.h

    def Solve(self,M1, M2, Tw01, Gamma, Eps):
        Init = 1.0
        T21 = Init
        OldTmp = T21
        Mbarsq = ((M1+M2)/2.0)**2.0
        Ft0 = self.FuncT0(Mbarsq, Gamma)
        Ff = self.FuncFf(Mbarsq, Gamma)
        Iter = 0
        while (True):
            Iter = Iter + 1
            Soln = -2.0 * Ff * ( np.log10(Tw01-1.0) - np.log10(Tw01 - T21))
            Soln = (M2**2.0 - M1**2.0) / Mbarsq + Soln
            Soln = 1.0 + 1.0 / Ft0 * Soln

            if (abs(Soln-OldTmp)<Eps):
                break
            if (Iter>100):
                Soln = (OldTmp+Soln)/2.0
                break
            
            OldTmp = Soln
            T21 = Soln
        


    def FuncT0(self, Mbarsq, Gamma):
        FuncT0 = (1.0+Gamma*Mbarsq)*(1.0+(Gamma-1.0)/2.0*Mbarsq)/(1.0-Mbarsq)
        return FuncT0

    def FuncFf(self, Mbarsq, Gamma):
        FuncFf = Gamma*Mbarsq*(1.0+(Gamma-1.0)/2.0*Mbarsq)/(1.0-Mbarsq)
        return FuncFf
