import numpy as np
from scipy import special

def calculate_N(x):
    output = 1 + special.erf(x/np.sqrt(2))
    return (0.5*output)

class Option:
#    all = []
    
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:

        assert t>=0, "time can not be negative"
        assert T>=0, "time can not be negative"
        assert r >= 0, "risk free rate can not be negative"

        self.S = S
        self.t = t
        self.K = K
        self.T = T
        self.sigma = sigma
        self.r = r

    @property
    def calculate_dPlus(self):
        dPlus = np.log(self.S/self.K) + (self.r + (0.5*self.sigma**2))*(self.T-self.t)
        dPlus /= (self.sigma*np.sqrt(self.T-self.t))
        return dPlus

    def calculate_dMinus(self):
        dMinus = np.log(self.S/self.K) - (self.r + (0.5*self.sigma**2))*(self.T-self.t)
        dMinus /= (self.sigma*np.sqrt(self.T-self.t))
        return dMinus

    def price(self): pass

    def delta(self): pass
    
    def get_delta(self):
        dPlus = self.calculate_dPlus()
        output = calculate_N(dPlus)
        return output 

    def __repr__(self) -> str:
        return f"C({self.t};{self.K},{self.T})"

class Call(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)

    def price(self):
        dPlus = self.calculate_dPlus()
        dMinus = self.calculate_dMinus()
        output = self.S*calculate_N(dPlus)
        output -= self.K*np.exp(-self.r*(self.T-self.t))*calculate_N(dMinus)
        return output

class Put(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)

    def price(self):
        dPlus = self.calculate_dPlus()
        dMinus = self.calculate_dMinus()
        output = self.K*np.exp(-self.r*(self.T-self.t))*calculate_N(-dMinus)
        output -= self.S*calculate_N(-dPlus)
        return output
