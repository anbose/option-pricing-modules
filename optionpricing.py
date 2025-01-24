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
        self.name = 'option'

    def calculate_dPlus(self):
        dPlus = np.log(self.S/self.K) + (self.r + (0.5*self.sigma**2))*(self.T-self.t)
        dPlus /= (self.sigma*np.sqrt(self.T-self.t))
        return dPlus

    def calculate_dMinus(self):
        dMinus = np.log(self.S/self.K) - (self.r + (0.5*self.sigma**2))*(self.T-self.t)
        dMinus /= (self.sigma*np.sqrt(self.T-self.t))
        return dMinus

    def delta(self): pass
    
    def get_delta(self):
        dPlus = self.calculate_dPlus()
        output = calculate_N(dPlus)
        return output

    def price(self):
        if self.name == 'Call':
            dPlus = self.calculate_dPlus()
            dMinus = self.calculate_dMinus()
            output = self.S*calculate_N(dPlus)
            output -= self.K*np.exp(-self.r*(self.T-self.t))*calculate_N(dMinus)
        elif self.name == 'Put':
            dPlus = self.calculate_dPlus()
            dMinus = self.calculate_dMinus()
            output = self.K*np.exp(-self.r*(self.T-self.t))*calculate_N(-dMinus)
            output -= self.S*calculate_N(-dPlus)
        else:
            raise ValueError('Option type not mentioned')
        
        return output
            

    def calculate_option_prices_for_parameters(self,parameter_name,parameter_range):
        option_prices = []
        for param in parameter_range:
            if parameter_name == 'volatility':
                self.sigma = param
                price = self.price()
            elif parameter_name == 'time':
                self.t = param
                price = self.price()
            elif parameter_name == 'underlying_price':
                self.S = param
                price = self.price()
            elif parameter_name == 'strike_price':
                self.K = param
                price = self.price()
            elif parameter_name == 'interest_rate':
                self.r = param
                price = self.price()
            elif parameter_name == 'maturity_time':
                self.T = param
                if self.t == self.T:
                    price = max(self.S-self.K,0) if self.name == 'Call' else max(self.K-self.S,0)
                else:
                    price = self.price()
            else:
                raise ValueError('Invalid parameter name. Possible parameter names are "volatility" : sigma, "time" : t, "underlying_price" : S, "strike_price" : K, "interest_rate" : r, "maturity_time" : T')

            option_prices.append(price)

        return option_prices

    def __repr__(self) -> str:
        return f"C({self.t};{self.K},{self.T})"

class Call(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)
        self.name = 'Call'

    #def price_vs_volatility(   
    #def plot_vs_time(self,ti,tf):
        

class Put(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)
        self.name = 'Put'

