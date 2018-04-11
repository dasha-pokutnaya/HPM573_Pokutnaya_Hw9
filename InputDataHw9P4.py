# simulation settings
POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 50    # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
ADD_BACKGROUND_MORT = True  # if background mortality should be added
DELTA_T = 1       # years

PSA_ON = False      # if probabilistic sensitivity analysis is on

# transition matrix
# Well Stroke Post-Stroke Dead
# What happens in the next year
TRANS_MATRIX = [
    [0.7375,  0.1575,    0.0,    0.105],   # Well
    [0,     0.0,    1,    0.0],   # Stroke
    [0,     0.20,      0.65,   0.15],   # Post-stroke
    ]

TREATMENT_RR = 1.05