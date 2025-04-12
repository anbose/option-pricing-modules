import streamlit as st
import optionpricing as opm
import numpy as np
import matplotlib.pyplot as plt

greek_dict = {
    "Delta" : "$\\Delta$",
    "Gamma" : "$\\Gamma$",
    "Vega" : "$V$",
    "Rho" : "$\\rho$",
}

param_dict = {
    "Underlying Price" : "underlying_price",
    "Strike Price" : "strike_price",
    "Time to Maturity (Years)" : "maturity_time",
    "Volatility" : "volatility",
    "Risk-free Interest Rate" : "interest_rate"
}

st.title('European Option Pricing and Greeks')
#st.divider()

## Sidebar (option pricing)

st.sidebar.header('Option Parameters',divider="rainbow")
S = st.sidebar.number_input("Underlying Asset Price (S)", value=100.0)
#t = st.sidebar.number_input("Time (current)", value=0.0)
K = st.sidebar.number_input("Strike Price (K)", value=110.0)
T = st.sidebar.number_input('Time to Maturity (T) in years', value=1.0)
sigma = st.sidebar.number_input("Volatility ($\\sigma$)", value=0.4)
r = st.sidebar.number_input('Risk-free Interest Rate (r)', value=0.3)

option_type = st.sidebar.selectbox("Option Type", ["Call", "Put"])

#price = 0
if st.sidebar.button("Calculate Price"):
    if option_type == "Call":
        price = opm.Call(S, 0, K, T, sigma, r).price()
    else:
        price = opm.Put(S, 0, K, T, sigma, r).price()
    st.sidebar.subheader(f"{option_type} Option Price: Euro {price:.2f}",divider='blue')

greek_name = st.sidebar.selectbox("Select Greek", ["Delta", "Gamma", "Vega", "Theta", "Rho"])
if st.sidebar.button('Calculate Greek Value'):
    if option_type == 'Call':
        greek_value = opm.Call(S, 0, K, T, sigma, r).calculate_greek(greek_name=greek_name)
    else:
        greek_value = opm.Put(S, 0, K, T, sigma, r).calculate_greek(greek_name=greek_name)
    st.sidebar.subheader(f'Option Greek Value: Euro {greek_value:.2f}',divider=True)



########################################################

st.header("Interactive Option Price Plot",divider='rainbow')

primary_parameter = st.selectbox("Select Primary Parameter to Plot", ["Underlying Price", "Strike Price", "Time to Maturity (Years)", "Volatility", 
                                                                      "Risk-free Interest Rate"])
col1, col2, col3 = st.columns(3)  # Create 3 columns

with col1:
    pg_i = st.number_input(f"Initial Value ({primary_parameter})", value=0.1,key=1)

with col2:
    pg_f = st.number_input(f"Final Value ({primary_parameter})", value=1.0,key=2)

#d_pg = 0.1 if pg_f >= 0.1 else 0.01

with col3:
    d_pg = st.number_input(f"Step Value ({primary_parameter})", value=0.1,key=3)

primary_data = np.arange(pg_i,pg_f,0.1)


# --- Create Streamlit sliders for secondary parameters ---

secondary_params = {
    'underlying_price' : {'label': 'Underlying Price', 'min': 0.1, 'max': 10000., 'value': 1000.0},
    'time': {'label': 'Time to Maturity (Years)', 'min': 0.1, 'max': 5.0, 'value': 1.0},
    'strike_price': {'label': 'Strike Price', 'min': 500., 'max': 10000., 'value': 1100.0},
    'volatility': {'label': 'Volatility', 'min': 0.2, 'max': 0.6, 'value': 0.3},
    'interest_rate': {'label': 'Risk-free Interest Rate', 'min': 0.02, 'max': 0.2, 'value': 0.03},
}

option_type = st.selectbox("Option Type", ["Call", "Put"],key=10)
if option_type == "Call":
    option_instance = opm.Call(100.,0.,1000.,10.,0.4,0.03) # Instantiate with some initial parameters
else:
    option_instance = opm.Put(100.,0.,1000.,10.,0.4,0.03)

st.subheader('Other Parameters',divider='rainbow')

slider_values = {}
for param_name, config in secondary_params.items():
    if config['label'] != primary_parameter:
        slider_values[param_name] = st.slider(
            config['label'],
            min_value=config['min'],
            max_value=config['max'],
            value=config['value'],
            step=0.01 if isinstance(config['min'], float) else 1
        )


# --- Call the modified plotting function with slider values ---
fig = option_instance.plot_price_with_streamlit_sliders(
    param_dict[primary_parameter],
    primary_data,
    slider_values
)

# --- Render the plot in Streamlit ---
st.pyplot(fig)
#st.plotly_chart(fig)

# --- Optional: Reset button ---
#if st.button("Reset Parameters"):
#    st.info("Reset functionality needs implementation using st.session_state for persistent changes.")
