
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

class pipe:
    def __init__(self,
                pipeLength,
                initialMach,
                wallTemp,
                gamma,
                T02_T01 = 1.05):
        self.L = pipeLength
        self.M1 = initialMach
        self.OrgM1 = initialMach
        self.finalMach = initialMach
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
        #self.findT02()
        self.T02_T01 = T02_T01
        self.T01 = 300

    def printResults(self):
        for i in range(len(self._M2)):
            print("{:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f}".format(self._M2[i],self._T0_T01[i],self._4fx_D[i],self._T_T1[i],self._P_P1[i],self._P0_P01[i]))

    def findT02(self):
        OrgM1 = self.M1
        M2 = self.M1 + self.h
        Prod = 1.0
        CommG = (self.gamma - 1.0)/2.0
        while (round(M2,5)<=1):
            Soln = self.Solve(self.M1,M2,self.Tw01,self.gamma,self.Eps)
            if Soln>1.0:
                Prod = Prod * Soln
                X = 2.0 * (np.log(self.Tw01-1.0)-np.log(self.Tw01-Prod))
                T = Prod * (1.0+CommG*OrgM1*OrgM1)/(1.0+CommG*M2*M2)
                P = OrgM1/M2*T**0.5
                P0 = Prod**0.5*OrgM1/M2*((1.0+CommG*M2*M2)/(1.0+CommG*OrgM1*OrgM1))**((self.gamma+1.0)/(2.0*(self.gamma-1.0)))
                #print("{:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f} , {:.4f}".format(M2,Prod,X,T,P,P0))
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

    def findM2(self):
        M2 = self.M1 + self.h
        Prod = 1.0
        # if(M2 > 1.0):
        #     print("There is no solution for the given input, \n The solution does not converge after M2 = {:.4f}".format(M2))
        # else:
        while True:
            if (M2 > 1.0):
                print("There is no solution for the given input, \n The solution does not converge after M2 = {:.4f}".format(M2))
                break
            else:
                Soln = self.Solve(self.M1,M2,self.Tw01,self.gamma,self.Eps)
                if Soln>1.0:
                    Prod = Prod * Soln
                    if (Prod-self.T02_T01 > 0.0):
                        Prod = Prod / Soln
                        self.Improv(self.M1,self.T02_T01,self.Tw01,self.gamma,self.Eps,self.h,Prod)
                        break
                    else:
                        self.M1 = M2
                        M2 = M2 + self.h
                else:
                    print("There is no solution for the given input, \n The solution does not converge after M2 = {:.4f}".format(M2))
                    break
        #print(self.finalMach)
        return self._M2[-1]
    
    # def constantTemp(self, M2):
    #     Mbarsq = ((self.M1 + M2) / 2.0)**2.0
    #     val = M2**2.0 - self.M1**2.0 - 2.0*(self.T02_T01 - 1.0) * (self.FuncT0(Mbarsq, self.gamma)/(self.T02_T01+1.0))+2.0*self.FuncFf(Mbarsq,self.gamma)/(2.0*self.Tw01-self.T02_T01-1.0)
    #     return val

    # def constantTemp(self):
    #     T02_T01 = -np.exp(np.log(self.Tw01 - 1.0) - 2.0 * self.f * self.L / self.D ) + self.Tw01
    #     T

    # def findM2fsolve(self):
    #     root = fsolve(self.constantTemp,0.4)
    #     print(root)#"rootsq = {:.6f} , root = {:.6f}".format(rootsq,root))
    #     #return root

    def Improv(self,M1,T02_T01,Tw01,gamma,Eps,h,Prod):
        NewM1 = M1
        Step = h/10.0
        NewM2 = M1 + Step
        while(Step > Eps):
            Soln = self.Solve(NewM1,NewM2,Tw01,gamma,Eps)
            Prod = Prod * Soln
            if (abs(Prod-T02_T01) < Eps):
                #print(NewM2)
                break
            else:
                if(Prod-T02_T01 > 0.0):
                    NewM2 = NewM2 - Step
                    Step = Step / 10.0
                    if Step<Eps:
                        #print(NewM2)
                        break
                    else:
                        NewM2 = NewM2 + Step
        CommG = (self.gamma - 1.0)/2.0
        X = 2.0 * (np.log(self.Tw01-1.0)-np.log(self.Tw01-Prod))
        T = Prod * (1.0+CommG*self.OrgM1*self.OrgM1)/(1.0+CommG*NewM2*NewM2)
        P = self.OrgM1/NewM2*T**0.5
        P0 = Prod**0.5*self.OrgM1/NewM2*((1.0+CommG*NewM2*NewM2)/(1.0+CommG*self.OrgM1*self.OrgM1))**((self.gamma+1.0)/(2.0*(self.gamma-1.0)))
        self._M2.append(NewM2)
        self._T0_T01.append(Prod)
        self._4fx_D.append(X)
        self._T_T1.append(T)
        self._P_P1.append(P)
        self._P0_P01.append(P0)
        print("M2 = {:.4f}".format(self._M2[-1]))
        print("T0/T01 = {:.4f}".format(self._T0_T01[-1]))
        print("4fx_D = {:.4f}".format(self._4fx_D[-1]))
        print("T_T1 = {:.4f}".format(self._T_T1[-1]))
        print("P_P1 = {:.4f}".format(self._P_P1[-1]))
        print("P0_P01 = {:.4f}".format(self._P0_P01[-1]))

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
            #define a relation between Tw01 and T21 for constant heat flux Tw01 = f(T21) or more generally Twn_T0n = f(T21,Tw1_T01)
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
#      a.findT02()
#     plt.plot(a.fx4_D(),a.M2(),label='Tw/T01= {:.1f}'.format(i/10.0))
# plt.ylabel('M2')
# plt.xlabel("$\\frac{4fx}{D}$")
# plt.legend()
# plt.grid()
# plt.show()

# for i in range(10,71,10):
#     print("Tw/T0 = ", i/10.0)
#     a = pipe(pipeLength=10,initialMach= 0.2,wallTemp = i/10.0,gamma = 1.4)
#     a.findT02()
#     plt.plot(a.fx4_D(),a.T_T1(),label='Tw/T01= {:.1f}'.format(i/10.0))
# plt.ylabel('Tw/T01')
# plt.xlabel("$\\frac{4fx}{D}$")
# plt.legend()
# plt.grid()
# plt.show()

# for i in range(30,41,10):
#     print("Tw/T0 = ", i/10.0)
#     a.findT02()
#     a = pipe(pipeLength=10,initialMach= 0.4,wallTemp = i/10.0,gamma = 1.4)
#     

# a = pipe(pipeLength=10,initialMach= 0.4,wallTemp = 4,gamma = 1.4, T02_T01 = 1.00575)
# a.findM2()

# a = pipe(pipeLength=10,initialMach= 0.4,wallTemp = 4,gamma = 1.4, T02_T01 = 1.00575)
# a.findM2fsolve()