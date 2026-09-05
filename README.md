# Drug Shelf-Life Simulator (Arrhenius Kinetics)

A Python tool that models chemical degradation kinetics and applies the Arrhenius equation to predict shelf life across storage temperatures (25°C/40°C/60°C), reflecting the same accelerated stability testing principles used industry-wide in pharmaceutical shelf-life estimation. Built to strengthen understanding of reaction kinetics while practicing NumPy and Matplotlib.

## What it does
- Simulates concentration/degradation vs. time for reaction orders 0, 1, or 2
- Applies the Arrhenius equation to compute rate constants at multiple temperatures
- Predicts shelf life (time to reach a defined degradation limit, e.g. 10%)
- Generates an Arrhenius plot (ln k vs 1/T) — the standard method for
  extrapolating from fast, hot accelerated data to slow, real-world
  room-temperature storage conditions
- Outputs a summary table and two plots (PNG)

## Tech used
Python, NumPy, Matplotlib

## How to run
Outputs `degradation_vs_time.png`, `arrhenius_plot.png`, and a printed shelf-life summary table.

## Example output
At 25°C, 40°C, and 60°C storage, the model predicts shelf life (time to 10% degradation) of roughly 724 days, 106 days, and 12 days respectively — showing the dramatic acceleration of degradation with temperature, consistent with real-world drug stability behavior.

