
# simulation settings
POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 50    # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
ADD_BACKGROUND_MORT = True  # if background mortality should be added
DELTA_T = 1       # years

PSA_ON = False      # if probabilistic sensitivity analysis is on

# transition matrix
TRANS_MATRIX = [
    [0.75,  0.15,    0.0,    0.1],   # Well
    [0,     0.0,    1,    0.0],   # Stroke
    [0,     0.25,      0.55,   0.2],   # Post-stroke
    ]

# treatment relative risk
TREATMENT_RR = 0.509
#TREATMENT_RR_CI = 0.365, 0.71  # lower 95% CI, upper 95% CI

