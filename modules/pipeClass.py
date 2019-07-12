#import pint
import numpy as np
import matplotlib.pyplot as plt
#u = pint.UnitRegistry()
#from scipy.integrate import solve_ivp

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
        self._M2 = [initialMach]
        self._T0_T01 = [1]
        self._4fx_D = [0]
        self._T_T1 = [1]
        self._P_P1 = [1]
        self._P0_P01 = [1]
        self.findT02()

    def findT02(self):
        OrM1 = self.M1
        M2 = self.M1 + self.h
        Prod = 1.0
        CommG = (self.gamma - 1.0)/2.0
        while (M2<=1):
            Soln = self.Solve(self.M1,M2,self.Tw01,self.gamma,self.Eps)
            if Soln>1.0:
                Prod = Prod * Soln
                X = 2.0 * (np.log(self.Tw01-1.0)-np.log(self.Tw01-Prod))
                T = Prod * (1.0+CommG*OrM1*OrM1)/(1.0+CommG*M2*M2)
                P = OrM1/M2*T**0.5
                P0 = Prod**0.5*OrM1/M2*((1.0+CommG*M2*M2)/(1.0+CommG*OrM1*OrM1))**((self.gamma+1.0)/(2.0*(self.gamma-1.0)))
                print("{:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f}".format(M2,Prod,X,T,P,P0))
                self._M2.append(M2)
                self._T0_T01.append(Prod)
                self._4fx_D.append(X)
                self._T_T1.append(T)
                self._P_P1.append(P)
                self._P0_P01.append(P0)
            else:
                print("solution is not converging after M2 = {:.4f}".format(M2))
                break
            self.M1 = M2
            M2 = M2 + self.h

    def M2(self):
        return self._M2
    
    def T0_T01(self):
        return self._T0_T01

    def fx4_D(self):
        return self._4fx_D

    def T_T1(self):
        return self._T_T1

    def P_P1(self):
        return self._P_P1
    
    def P0_P01(self):
        return self._P0_P01

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
            Soln = -2.0 * Ff * ( np.log(Tw01-1.0) - np.log(Tw01 - T21))
            Soln = (M2**2.0 - M1**2.0) / Mbarsq + Soln
            Soln = 1.0 + 1.0 / Ft0 * Soln

            if (abs(Soln-OldTmp)<Eps):
                break
            if (Iter>100):
                Soln = (OldTmp+Soln)/2.0
                break
            if (Soln>Tw01):
                Soln = abs(Tw01-Soln)
            OldTmp = Soln
            T21 = Soln
        return Soln
        


    def FuncT0(self, Mbarsq, Gamma):
        FuncT0 = (1.0+Gamma*Mbarsq)*(1.0+(Gamma-1.0)/2.0*Mbarsq)/(1.0-Mbarsq)
        return FuncT0

    def FuncFf(self, Mbarsq, Gamma):
        FuncFf = Gamma*Mbarsq*(1.0+(Gamma-1.0)/2.0*Mbarsq)/(1.0-Mbarsq)
        return FuncFf

# for i in range(10,100, 5):
#     print("Tw/T0 = ", i/10.0)
#     a = pipe(pipeLength=10,initialMach= 0.5,wallTemp = i/10.0,gamma = 1.4)
#     plt.plot(a.fx4_D(),a.M2(),label='Tw/T01= {:.1f}'.format(i/10.0))
# plt.ylabel('M2')
# plt.xlabel("$\\frac{4fx}{D}$")
# plt.legend()
# plt.grid()
# plt.show()

for i in range(30,31):
    print("Tw/T0 = ", i/10.0)
    a = pipe(pipeLength=10,initialMach= 0.4,wallTemp = i/10.0,gamma = 1.4)
    plt.plot(a.fx4_D(),a.T_T1(),label='Tw/T01= {:.1f}'.format(i/10.0))
plt.ylabel('Tw/T01')
plt.xlabel("$\\frac{4fx}{D}$")
plt.legend()
plt.grid()
plt.show()
