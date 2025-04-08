# European Option Pricing and Greeks Web Application

This repository contains a Python web application built using Streamlit to calculate European option prices based on the Black-Scholes formula and visualize various option Greeks interactively.

## Overview

This application allows users to input parameters for European call and put options (such as underlying asset price, strike price, time to maturity, risk-free interest rate, and volatility) and instantly see the calculated option price. Additionally, it provides interactive visualizations of the option price against multiple parameters, with sliders to adjust other parameters in real-time.

## Key Features

* **Black-Scholes Pricing:** Calculates the price of European call and put options using the Black-Scholes formula.
* **Interactive Input:** User-friendly interface to input option parameters.
* **Option Greek Visualization:** Generates interactive plots of option price against a selected parameter.
* **Real-time Parameter Adjustment:** Utilizes Streamlit sliders to dynamically change secondary parameters and observe the impact on the Greek plots.
* **Clear Output:** Displays the calculated option price and renders interactive plots directly in the web browser.

## Prerequisites

* Python 3.6 or higher
* Git (for cloning the repository)
* pip (Python package installer)

## Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/anbose/option-pricing-modules.git
    ```
    
2.  **Navigate to the project directory:**
    ```bash
    cd option-pricing-module
    ```
    
3.  **Install dependencies:**
    It's recommended to create a virtual environment to manage dependencies. You can do this using `venv` or `conda`.
    

## How to Run

1.  Ensure you are in the project directory (where `app.py` is located) and your virtual environment is activated (if you created one).
2.  Run the Streamlit application using the following command:
    ```bash
    streamlit run app.py
    ```
3.  The application will automatically open in your default web browser. If it doesn't, you can usually find the URL in your terminal output (it's often `http://localhost:xxxx`).

## Usage

1.  **Option Price Calculation:**
    * In the left sidebar, you will find input fields for the option parameters: Underlying Asset Price (S), Strike Price (K), Time to Maturity (T), Risk-Free Interest Rate (r), and Volatility (σ).
    * Select the "Option Type" (Call or Put) using the dropdown menu.
    * Click the "Calculate Price" button. The calculated option price will be displayed below the button.
    * You will also find a dropdown menu below, to calculate any of the Greeks. Select one.
    * Click the "Calculate Greek Value" button. The calculated option greek will be displayed below the button.

2.  **Interactive Visualization:**
    * In the main section of the application, you will find a "Interactive Option Price Plot" header.
    * Use the "Select Primary parameter" dropdown to choose the primary parameter. The option price will be plotted against this parameter.
    * Set the range of the primary parameter using the input fields for "Initial Value (parameter name)", "Final Value (parameter name)" and "Step Value (parameter name)".
    * Select the "Option Type" using the dropdown menu. 
    * **Interactive Sliders:** Above the plot, you will find sliders for the secondary option parameters (e.g., Strike Price, Volatility, Interest Rate). You can adjust these sliders in real-time to observe how the option price changes with different parameter values.

## Contributing

If you'd like to contribute to this project, please feel free to:

1.  Fork the repository.
2.  Create a new branch for your feature or bug fix.
3.  Make your changes and commit them.
4.  Push your changes to your fork.
5.  Submit a pull request.

Please ensure your code follows the project's coding style and includes appropriate tests.

---

Enjoy using the European Option Pricing Web Application!
