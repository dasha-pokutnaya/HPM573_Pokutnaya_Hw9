from enum import Enum
import numpy as np
import scipy.stats as stat
import math as math
import InputDataHw9P3 as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.ProbDistParEst as Est


class HealthStats(Enum):
    """ health states of patients with HIV """
    Well  = 0
    Stroke = 1
    PostStroke = 2
    Dead = 3

class Therapies(Enum):
    #coagulation
    coagulation = 1


class _Parameters:

    def __init__(self, therapy, coagulation):

        # simulation time step
        self._delta_t = Data.DELTA_T

        self._therapy = therapy

        # initial health state
        self._initialHealthState = HealthStats.Well

        # transition probability matrix of the selected therapy
        self._prob_matrix = []

    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t

    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]


class ParametersFixed(_Parameters):
    def __init__(self, therapy):

        # initialize the base class
        _Parameters.__init__(self, therapy)

        # calculate transition probabilities between hiv states
        self._prob_matrix = calculate_prob_matrix()

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == therapy:
            # treatment relative risk
            self._treatmentRR = Data.TREATMENT_RR
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = calculate_prob_matrix_combo(
                matrix_mono=self._prob_matrix, combo_rr=Data.TREATMENT_RR)


class ParametersProbabilistic(_Parameters):
    def __init__(self, seed, therapy):

        # initializing the base class
        _Parameters.__init__(self, therapy)

        self._rng = Random.RNG(seed)    # random number generator to sample from parameter distributions
        self._lnRelativeRiskRVG = None  # random variate generator for the treatment relative risk
        self._annualStateCostRVG = []       # list of random variate generators for the annual cost of states
        self._annualStateUtilityRVG = []    # list of random variate generators for the annual utility of states

        # HIV transition probabilities
        #j = 0
        #for prob in Data.TRANS_MATRIX:
        #    self._hivProbMatrixRVG.append(Random.Dirichlet(prob[j:]))
        #    j += 1

        # treatment relative risk
        # find the mean and st_dev of the normal distribution assumed for ln(RR)
        sample_mean_lnRR = math.log(Data.TREATMENT_RR)
        sample_std_lnRR = (Data.TREATMENT_RR_CI[1]-Data.TREATMENT_RR_CI[0])/(2*stat.norm.ppf(1-0.05/2))
        self._lnRelativeRiskRVG = Random.Normal(mean=sample_mean_lnRR, st_dev=sample_std_lnRR)

        # annual state cost
        for cost in Data.ANNUAL_STATE_COST:
            # find shape and scale of the assumed gamma distribution
            shape, scale = Est.get_gamma_parameters(mean=cost, st_dev=cost/ 4)
            # append the distribution
            self._annualStateCostRVG.append(Random.Gamma(shape, scale))

        # resample parameters
        self.__resample()

    def __resample(self):

        # calculate transition probabilities
        # create an empty matrix populated with zeroes
        self._prob_matrix = []
        for s in HealthStats:
            self._prob_matrix.append([0] * len(HealthStats))

        # for all health states
        for s in HealthStats:
            # if the current state is death
            if s in [HealthStats.HIV_DEATH, HealthStats.BACKGROUND_DEATH]:
                # the probability of staying in this state is 1
                self._prob_matrix[s.value][s.value] = 1
            else:
                # sample from the dirichlet distribution to find the transition probabilities between hiv states
                dist = self._hivProbMatrixRVG[s.value]
                sample = dist.sample(self._rng)
                for j in range(len(sample)):
                    self._prob_matrix[s.value][s.value+j] = sample[j]

        # sample from gamma distributions that are assumed for annual state costs
        self._annualStateCosts = []
        for dist in self._annualStateCostRVG:
            self._annualStateCosts.append(dist.sample(self._rng))

        # sample from beta distributions that are assumed for annual state utilities
        self._annualStateUtilities = []
        for dist in self._annualStateUtilityRVG:
            self._annualStateUtilities.append(dist.sample(self._rng))


def calculate_prob_matrix():
    """ :returns transition probability matrix for stroke states"""

    # create an empty matrix populated with zeroes
    prob_matrix = []
    for s in HealthStats:
        prob_matrix.append([0] * len(HealthStats)) #there are 4 health states

    # for all health states
    for s in HealthStats:
        # if the current state is death
        if s in [HealthStats.Dead]:
            # the probability of staying in this state is 1
            prob_matrix[s.value][s.value] = 1


    return prob_matrix

# we counted how many well individuals started and how many people
# ended in a different state


def calculate_prob_matrix_combo(matrix_mono, combo_rr):
    """
    :param matrix_mono: (list of lists) transition probability matrix under mono therapy
    :param combo_rr: relative risk of the combination treatment
    :returns (list of lists) transition probability matrix under combination therapy """

    # create an empty list of lists
    matrix_combo = []
    for l in matrix_mono:
        matrix_combo.append([0] * len(l))

    # populate the combo matrix
    # first non-diagonal elements
    for s in HealthStats:
        for next_s in range(s.value + 1, len(HealthStats)):
            matrix_combo[s.value][next_s] = combo_rr * matrix_mono[s.value][next_s]

    # diagonal elements are calculated to make sure the sum of each row is 1
    for s in HealthStats:
        if s not in [HealthStats.Dead, HealthStats.BACKGROUND_DEATH]:
            matrix_combo[s.value][s.value] = 1 - sum(matrix_combo[s.value][s.value + 1:])

    return matrix_combo
