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
    "Time to Maturity" : "maturity_time",
    "Volatility" : "volatility",
    "Risk-free Interest Rate" : "interest_rate"
}

st.title('European Option Pricing and Greeks')

st.sidebar.header('Option Parameters')
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
    st.subheader(f"{option_type} Option Price: Euro {price:.2f}")

greek_name = st.sidebar.selectbox("Select Greek", ["Delta", "Gamma", "Vega", "Theta", "Rho"])
if st.sidebar.button('Calculate Greek Value'):
    if option_type == 'Call':
        greek_value = opm.Call(S, 0, K, T, sigma, r).calculate_greek(greek_name=greek_name)
    else:
        greek_value = opm.Put(S, 0, K, T, sigma, r).calculate_greek(greek_name=greek_name)
    st.subheader(f'Option Greek {greek_name} ({greek_dict[greek_name]}) Value: Euro {greek_value:.2f}')



########################################################

st.header("Interactive Option Price Plot")

primary_parameter = st.selectbox("Select Primary Parameter to Plot", ["Underlying Price", "Strike Price", "Time to Maturity", "Volatility", 
                                                                      "Risk-free Interest Rate"])
col1, col2, col3 = st.columns(3)  # Create 3 columns

with col1:
    pg_i = st.number_input(f"Initial Value ({primary_parameter})", value=0.1,key=1)

with col2:
    pg_f = st.number_input(f"Final Value ({primary_parameter})", value=1.0,key=2)

with col3:
    d_pg = st.number_input(f"Step Value ({primary_parameter})", value=0.1,key=3)

primary_data = np.arange(pg_i,pg_f,d_pg)

#sec_greek1 = st.selectbox("Select First Secondary Greek to Plot", ["Delta", "Gamma", "Vega", "Theta", "Rho"])
#col3, col4, col5 = st.columns(3)  # Create 3 columns

#with col3:
#    sg1_i = st.number_input(f"Initial Value ({sec_greek1})", value=0.1,key=4)

#with col4:
#    sg1_f = st.number_input(f"Final Value ({sec_greek1})", value=1.0,key=5)

#with col5:
#    d_sg1 = st.number_input(f"Step Value ({sec_greek1})", value=0.1,key=6)

#sg1_data = np.arange(sg1_i,sg1_f,d_sg1)

#sec_greek2 = st.selectbox("Select Second Secondary Greek to Plot", ["Delta", "Gamma", "Vega", "Theta", "Rho"])
#col7, col8, col9 = st.columns(3)  # Create 3 columns

#with col7:
#    sg2_i = st.number_input(f"Minimum ({sec_greek2})", value=0.1,key=7)

#with col8:
#    sg2_f = st.number_input(f"Maximum ({sec_greek2})", value=1.0,key=8)

#with col9:
#    d_sg2 = st.number_input(f"Step ({sec_greek2})", value=0.1,key=9)

#sg2_data = np.arange(sg2_i,sg2_f,d_sg2)

# --- Create Streamlit sliders for secondary parameters ---
secondary_params = {
    #'time': {'label': 'Time to Maturity (Years)', 'min': 0.1, 'max': 5.0, 'value': 1.0},
    'strike_price': {'label': 'Strike Price', 'min': 500., 'max': 1500., 'value': 1100.0},
    'volatility': {'label': 'Volatility ($\\sigma$)', 'min': 0.2, 'max': 0.6, 'value': 0.3},
    'interest_rate': {'label': 'Interest Rate (r)', 'min': 0.02, 'max': 0.2, 'value': 0.03},
}


slider_values = {}
for param_name, config in secondary_params.items():
    slider_values[param_name] = st.slider(
        config['label'],
        min_value=config['min'],
        max_value=config['max'],
        value=config['value'],
        step=0.01 if isinstance(config['min'], float) else 1
    )

option_type = st.selectbox("Option Type", ["Call", "Put"],key=10)
if option_type == "Call":
    option_instance = opm.Call(100.,0.,1000.,10.,0.4,0.03) # Instantiate with some initial parameters
else:
    option_instance = opm.Put(100.,0.,1000.,10.,0.4,0.03)

# --- Call the modified plotting function with slider values ---
fig = option_instance.plot_price_with_streamlit_sliders(
    param_dict[primary_parameter],
    primary_data,
    slider_values
)

# --- Render the plot in Streamlit ---
st.pyplot(fig)

# --- Optional: Add a reset button ---
if st.button("Reset Parameters"):
    # Implement logic to reset slider values (using st.session_state might be needed for more complex resets)
    st.info("Reset functionality needs implementation using st.session_state for persistent changes.")
