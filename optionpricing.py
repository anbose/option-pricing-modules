import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import copy
from matplotlib.widgets import Button, Slider
#plt.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Computer Modern']


def calculate_N(x):
    output = 1 + special.erf(x/np.sqrt(2))
    return (0.5*output)

class Option:
    
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:

        assert t >= 0, "current time (t) can not be negative"
        assert T > t, "maturity time (T) must be greater than current time (t)"
        assert r >= 0, "risk free rate (r) can not be negative"
        assert sigma >= 0, "Volatility (sigma) can not be negative"

        self.S = S
        self.t = t
        self.K = K
        self.T = T
        self.sigma = sigma
        self.r = r
        self.name = 'option'

    def calculate_dPlus(self):
        log_term = np.log(self.S/self.K)
        time_diff = self.T - self.t
        dPlus = log_term + (self.r + (0.5*self.sigma**2))*time_diff
        dPlus /= (self.sigma*np.sqrt(time_diff))
        return dPlus

    def calculate_dMinus(self):
        log_term = np.log(self.S/self.K)
        time_diff = self.T - self.t
        dMinus = log_term + (self.r - (0.5*self.sigma**2))*time_diff
        dMinus /= (self.sigma*np.sqrt(time_diff))
        return dMinus

    def Delta(self):
        if self.name == 'Call':
            return calculate_N(self.calculate_dPlus())
        elif self.name == 'Put':
            return -calculate_N(-self.calculate_dPlus())
        else:
            raise ValueError('Option type not mentioned')

    def Gamma(self):
        dPlus = self.calculate_dPlus()
        Nprime = np.exp(-0.5*(dPlus**2)) / (np.sqrt(2*np.pi))
        time_diff = self.T - self.t
        return (Nprime / (self.S * self.sigma * np.sqrt(time_diff)))

    def Vega(self):
        dPlus = self.calculate_dPlus()
        Nprime = np.exp(-0.5*(dPlus**2)) / (np.sqrt(2*np.pi))
        time_diff = self.T - self.t
        output = self.S * Nprime * np.sqrt(time_diff)
        return output

    def Theta(self):
        dPlus = self.calculate_dPlus()
        dMinus = self.calculate_dMinus()
        Nprime = np.exp(-0.5*(dPlus**2)) / (np.sqrt(2*np.pi))
        time_diff = self.T - self.t
        output_part1 = - (self.S * Nprime * self.sigma) / (2 * np.sqrt(time_diff))
        if self.name == 'Call':
            output_part2 = self.r * self.K * calculate_N(dMinus) * np.exp(-self.r * time_diff)
            output = output_part1 - output_part2
        elif self.name == 'Put':
            output_part2 = self.r * self.K * calculate_N(-dMinus) * np.exp(-self.r * time_diff)
            output = output_part1 + output_part2
        else:
            raise ValueError('Option type not mentioned')

        return output
   
    def Rho(self):
        dMinus = self.calculate_dMinus()
        time_diff = self.T - self.t
        if self.name == 'Call':
            output = self.K * time_diff * np.exp(-self.r * time_diff) * calculate_N(dMinus)
        elif self.name == 'Put':
            output = - self.K * time_diff * np.exp(-self.r * time_diff) * calculate_N(-dMinus)
        else:
            raise ValueError('Option type not mentioned')

        return output

    greek_dict = {
        'Delta' : lambda option : option.Delta(),
        'Gamma' : lambda option : option.Gamma(),
        'Vega' : lambda option : option.Vega(),
        'Theta' : lambda option : option.Theta(),
        'Rho' : lambda option : option.Rho(),
    }

    def calculate_greek(self,greek_name):
        if greek_name in self.greek_dict:
            return (self.greek_dict[greek_name](self))
        else:
            raise ValueError('incorrect option greek. Possible greeks are "Delta", "Gamma", "Vega", "Theta" and "Rho"')
    
    def price(self):
        dPlus = self.calculate_dPlus()
        dMinus = self.calculate_dMinus()
        time_diff = self.T - self.t

        if self.name == 'Call':
            output = self.S * calculate_N(dPlus)
            output -= self.K * np.exp(-self.r * time_diff) * calculate_N(dMinus)
        elif self.name == 'Put':
            output = self.K * np.exp(-self.r * time_diff) * calculate_N(-dMinus)
            output -= self.S * calculate_N(-dPlus)
        else:
            raise ValueError('Option type not mentioned')
        
        return output
            
    def calculate_option_prices_for_parameters(self,parameter_name,parameter_range):
        option_prices = []
        
        option_copy = copy.deepcopy(self)

        param_setter = {
            'underlying_price': lambda option,val : setattr(option, 'S', val),
            'time': lambda option,val : setattr(option, 't', val),
            'strike_price': lambda option,val : setattr(option, 'K', val),
            'maturity_time': lambda option,val : setattr(option, 'T', val),
            'volatility': lambda option,val : setattr(option, 'sigma', val),
            'interest_rate': lambda option,val : setattr(option, 'r', val)
        }
        
        if parameter_name in param_setter:
            for param in parameter_range:
                param_setter[parameter_name](option_copy,param)
                if parameter_name == 'maturity_time' and option_copy.t == option_copy.T:
                    price = max(option_copy.S - option_copy.K,0) if option_copy.name == 'Call' else max(option_copy.K-option_copy.S,0)
                else:
                    price = option_copy.price()

                option_prices.append(price)
        else:
            raise ValueError('Invalid parameter name. Possible parameter names are  "underlying_price" : S, "time" : t, "strike_price" : K, "maturity_time" : T, "volatility" : sigma, "interest_rate" : r') 

        return option_prices

    def plot_price(self,parameter_name,parameter_range):
        ydata = self.calculate_option_prices_for_parameters(parameter_name,parameter_range)
        xdata = parameter_range
        parameter_label = ' '.join(parameter_name.split('_'))

        fig,ax = plt.subplots(1,1,figsize=(10,6))
        ax.set_title('price vs ' + parameter_label,fontsize=16)
        ax.plot(xdata,ydata,'o',markersize=3)
        ax.set_xlabel(parameter_label,fontsize=16)
        ax.set_ylabel('Option Price',fontsize=16)
        ax.grid(True)
        plt.tight_layout()
        plt.show()

    def plot_price_with_silders(self,*args):
        primary_parameter = args[0]
        primary_data = args[1]
        secondary_params = {}
        for param_id in range(2,len(args),2):
            param_name = args[param_id]
            param_data = args[param_id+1]
            secondary_params[param_name] = param_data

        param_setter = {
            'underlying_price': lambda option,val : setattr(option, 'S', val),
            'time': lambda option,val : setattr(option, 't', val),
            'strike_price': lambda option,val : setattr(option, 'K', val),
            'maturity_time': lambda option,val : setattr(option, 'T', val),
            'volatility': lambda option,val : setattr(option, 'sigma', val),
            'interest_rate': lambda option,val : setattr(option, 'r', val)
        }
        
        fig,ax = plt.subplots(figsize=(10,6))
        plt.subplots_adjust(left=0.25,bottom=0.25)
        price_data = self.calculate_option_prices_for_parameters(primary_parameter,primary_data)
        line, = ax.plot(primary_data,price_data,'o',markersize=3)
        ax.set_xlabel(' '.join(primary_parameter.split('_')),fontsize=16)
        ax.set_ylabel('Option Price',fontsize=16)
        ax.grid(True)

        slider_params = list(secondary_params.keys())
        slider_dict = {}
        for slider_id in range(len(secondary_params)):
            slider_ax = fig.add_axes([0.25,0.12-slider_id*0.05,0.65,0.03])      # [left, bottom, width, height]
            slider_label = ' '.join(slider_params[slider_id].split('_'))
            slider_range = secondary_params[slider_params[slider_id]]
            slider_dict[slider_params[slider_id]] = Slider(ax=slider_ax,
                                                            label=slider_label,valmin=slider_range.min(),
                                                                valmax=slider_range.max(),valinit=slider_range[0])
        current_params = {}
        def update(val):
            option_copy = copy.deepcopy(self)
            for params in slider_dict:
                param_setter[params](option_copy,slider_dict[params].val)

            prices = option_copy.calculate_option_prices_for_parameters(primary_parameter,primary_data)

            line.set_ydata(prices)
            ax.relim()
            ax.autoscale_view()
            fig.canvas.draw_idle()

        for params in slider_dict:
            slider_dict[params].on_changed(update)

        reset_ax = fig.add_axes([0.1,0.25,0.05,0.03])
        reset_button = Button(reset_ax,'Reset',color='0.85',hovercolor='0.95',useblit=True)

        def reset(event):
            for params in slider_dict:
                slider_dict[params].reset()
        reset_button.on_clicked(reset)

        plt.show()

    def plot_greek(self,greek_name,param_name,param_data):
        greeks = []

        option_copy = copy.deepcopy(self)
        param_setter = {
            'underlying_price': lambda option,val : setattr(option, 'S', val),
            'time': lambda option,val : setattr(option, 't', val),
            'strike_price': lambda option,val : setattr(option, 'K', val),
            'maturity_time': lambda option,val : setattr(option, 'T', val),
            'volatility': lambda option,val : setattr(option, 'sigma', val),
            'interest_rate': lambda option,val : setattr(option, 'r', val)
        }

        if param_name in param_setter:
            for param in param_data:
                param_setter[param_name](option_copy,param)
                greek_value = option_copy.calculate_greek(greek_name)
                greeks.append(greek_value)
        else:
            raise ValueError('Invalid parameter name. Possible parameter names are  "underlying_price" : S, "time" : t, "strike_price" : K, "maturity_time" : T, "volatility" : sigma, "interest_rate" : r')
        
        fig,ax = plt.subplots(figsize=(10,6))
        plt.plot(param_data,greeks,'o',markersize=3)
 
        param_label = ' '.join(param_name.split('_'))
        ax.set_title(f'Option Greek {greek_name} vs {param_label}')
        ax.set_xlabel(param_label)
        ax.set_ylabel(greek_name)
        ax.grid(True)
        plt.tight_layout()
        plt.show()

    def implied_volatility(self,p0,v0):
        option_copy = copy.deepcopy(self)
        
        assert v0 > 0, "volatility can not be negative"
        
        while True:
            setattr(option_copy,'sigma',v0)
            fx = option_copy.price() - p0
            v_next = option_copy.sigma - (fx/option_copy.Vega())
            if (abs(v_next - v0) < 1e-6):
                break
            v0 = v_next
        return v0

    def __repr__(self) -> str:
        return f"C({self.t};{self.K},{self.T})"

class Call(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)
        self.name = 'Call'

class Put(Option):
    def __init__(self,S:float,t:float,K:float,T:float,sigma:float,r:float) -> None:
        super().__init__(S,t,K,T,sigma,r)
        self.name = 'Put'

