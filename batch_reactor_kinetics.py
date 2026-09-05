"""
BATCH REACTOR KINETICS SIMULATOR
--------------------------------
What this does, in one line:
    It predicts how fast a chemical (e.g. a drug in a tablet) degrades
    over time, at different storage temperatures, and estimates shelf life.

Why this matters for pharma (Dr. Reddy's):
    Real pharma companies run "stability studies" where they store a drug
    at a few fixed temperatures (commonly 25C, 40C, 60C - this mirrors the
    ICH guideline conditions used industry-wide) and measure how much of
    the drug has degraded after some time. They then use the Arrhenius
    equation (the same one used here) to predict how the drug will behave
    at normal storage temperature over its full shelf life, without having
    to wait years to find out.

Author: Pari
"""

import numpy as np
import matplotlib.pyplot as plt

# ======================================================================
# STEP 1: INPUTS YOU CAN CHANGE
# ======================================================================
C0 = 100.0            # starting amount of drug, in % potency (100% = fresh)
order = 1              # reaction order (1 = first order; most drug
                       # degradation reactions are modeled as first order)

A_factor = 4.83e13     # Arrhenius pre-exponential factor, units: 1/day
Ea = 100000            # activation energy of the degradation reaction, J/mol
R_gas = 8.314          # universal gas constant, J/(mol*K) -- this is fixed,
                       # it is a physical constant, not something to change

temperatures_C = [25, 40, 60]   # storage temperatures to compare, deg C
                                 # (25C = long-term/room storage,
                                 #  40C = accelerated condition,
                                 #  60C = stress condition -- these match
                                 #  real ICH stability-testing conditions)

t_max = 800            # how many days to simulate
n_points = 400          # how many points to calculate along that time span
shelf_life_limit = 0.10  # drug is considered expired once 10% has degraded
                          # (a common real-world regulatory shelf-life rule)


# ======================================================================
# STEP 2: ARRHENIUS EQUATION -- how temperature affects reaction speed
# ======================================================================
def arrhenius_k(A_factor, Ea, T_kelvin):
    """
    k = A * exp( -Ea / (R * T) )

    In plain English: as temperature (T) goes up, the negative exponent
    becomes less negative, so exp(...) gets bigger, so k (the rate
    constant) gets bigger. Bigger k = faster degradation.
    """
    return A_factor * np.exp(-Ea / (R_gas * T_kelvin))


# ======================================================================
# STEP 3: INTEGRATED RATE LAW -- how concentration changes with time
# ======================================================================
def concentration_profile(C0, k, order, t):
    """
    Starting from the rate law  -dC/dt = k * C^order,
    solving that differential equation (integrating it) gives a direct
    formula for concentration C at any time t. That's what's below --
    we don't need to solve the equation live, just plug numbers into
    the already-solved formula.
    """
    if order == 0:
        C = C0 - k * t
        C = np.clip(C, 0, None)     # concentration can never go below 0
    elif order == 1:
        C = C0 * np.exp(-k * t)
    elif order == 2:
        C = C0 / (1 + k * C0 * t)
    else:
        raise ValueError("order must be 0, 1, or 2")
    return C


# ======================================================================
# STEP 4: RUN THE SIMULATION FOR EACH TEMPERATURE
# ======================================================================
t = np.linspace(0, t_max, n_points)   # an array of time values: 0, ..., t_max

plt.figure(figsize=(8, 6))

results_summary = []   # we'll store (temperature, k, shelf_life_days) here

for T_C in temperatures_C:
    T_K = T_C + 273.15                       # convert Celsius to Kelvin
    k = arrhenius_k(A_factor, Ea, T_K)       # rate constant at this temp
    C = concentration_profile(C0, k, order, t)
    fraction_degraded = 1 - (C / C0)         # 0 = nothing degraded, 1 = all gone

    plt.plot(t, fraction_degraded * 100,
             label=f"{T_C} deg C  (k = {k:.2e} /day)")

    # find the first day where degradation crosses the shelf-life limit
    over_limit = fraction_degraded >= shelf_life_limit
    if np.any(over_limit):
        shelf_life_days = t[np.argmax(over_limit)]
    else:
        shelf_life_days = None  # never reached within t_max

    results_summary.append((T_C, k, shelf_life_days))

plt.axhline(shelf_life_limit * 100, color="gray", linestyle="--",
            label=f"{int(shelf_life_limit*100)}% degradation limit")
plt.xlabel("Time (days)")
plt.ylabel("Degradation (%)")
plt.title("Drug Degradation vs Time at Different Storage Temperatures")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("degradation_vs_time.png", dpi=150)
plt.close()


# ======================================================================
# STEP 5: PRINT A SUMMARY TABLE
# ======================================================================
print("\nSHELF-LIFE SUMMARY (time to reach {}% degradation)".format(
    int(shelf_life_limit * 100)))
print(f"{'Temp (C)':>10} {'k (1/day)':>15} {'Shelf life (days)':>20}")
for T_C, k, shelf_life_days in results_summary:
    life_str = f"{shelf_life_days:.1f}" if shelf_life_days is not None else "not reached"
    print(f"{T_C:>10} {k:>15.3e} {life_str:>20}")


# ======================================================================
# STEP 6: ARRHENIUS PLOT (ln k vs 1/T)
#          This is exactly the plot pharma stability scientists use to
#          extrapolate from fast, hot, accelerated data down to slow,
#          cool, real-world storage conditions.
# ======================================================================
T_kelvin_arr = np.array([T_C + 273.15 for T_C in temperatures_C])
k_arr = np.array([arrhenius_k(A_factor, Ea, T) for T in T_kelvin_arr])

plt.figure(figsize=(6, 5))
plt.plot(1 / T_kelvin_arr, np.log(k_arr), "o-", color="darkred")
plt.xlabel("1 / T  (1/K)")
plt.ylabel("ln(k)")
plt.title("Arrhenius Plot")
plt.grid(True)
plt.tight_layout()
plt.savefig("arrhenius_plot.png", dpi=150)
plt.close()

print("\nSaved plots: degradation_vs_time.png and arrhenius_plot.png")
