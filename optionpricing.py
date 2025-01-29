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

    def delta(self):
        if self.name == 'Call':
            return calculate_N(self.calculate_dPlus())
        elif self.name == 'Put':
            return -calculate_N(-self.calculate_dPlus())
        else:
            raise ValueError('Option type not mentioned')
    
    def price(self):
        dPlus = self.calculate_dPlus()
        dMinus = self.calculate_dMinus()
        time_diff = self.T - self.t

        if self.name == 'Call':
            output = self.S * calculate_N(dPlus)
            output -= self.K * np.exp(-self.r * time_diff) * calculate_N(dMinus)
        elif self.name == 'Put':
            output = self.K * np.exp(-self.r* time_diff) * calculate_N(-dMinus)
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

    def plot_greek(self,parameter_name,parameter_range):
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

    def plot_with_silders(self,*args):
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

