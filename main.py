#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Babin.VM
#
# Created:     02.09.2018
# Copyright:   (c) Babin.VM 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import calendar
import copy
import datetime
import math
import random
import scipy
import os

import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as AA
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.optimize import minimize, rosen, rosen_der
from matplotlib.ticker import FuncFormatter
from itertools import permutations

data_dir = os.path.join(os.getcwd(), 'chatylka_data')
from chatylka_data.chatylka_constants import *

from src.common import *

from src.schedule import CreateSchedule
from src.report import Report
from src.well import GetWellsFromPreparedTxtFiles
from src.profiles import CalculateProfiles2

class OptimalSolution:
    def __init__(self, method, wells, strategy, Qmax):
        if len(wells) != len(strategy):
            raise ValueError('Wells count and strategy lists length are not equal!')

        self.method = method
        self.Qmax = Qmax
        self.optimal_strategy = strategy
        self.wells = wells

        self.Einj = sum(strategy)

        self.Q_wells = []
        self.I_wells = []
        self.P_wells = []

        for i, w in enumerate(wells):
            Q = float(w.Q_E_interpolate(self.optimal_strategy[i]))

            Q_closest = min(w.Q_I, key=lambda y: abs(y - Q))
            closest_index = (w.Q_I).index(Q_closest)

            if abs(closest_index - len(w.I)) <= 2:
                print('WARNING: for well ', w.name, ' injected volume = INF., but calculated as finite volume.')
                print('closest index = ', closest_index, ' out of ', len(w.I))

            I = (w.I)[closest_index]
            P = (w.P_I)[closest_index]

            self.Q_wells.append(round(Q, 3))
            self.I_wells.append(round(I, 3))
            self.P_wells.append(round(P, 3))

        # check 2
        Qmax_check = sum(self.Q_wells)
        if abs(Qmax_check - self.Qmax) / abs(self.Qmax) > 0.01:
            raise ValueError('inconsistent data in Solution, Qmax (external) = %f; Qmax (internal) = %f',
                             Qmax, Qmax_check)
        # end check 2

    def PrintInfo(self):
        print(str(self.method) + ' optimization result:')
        print('Q max = ', self.Qmax)
        print('Optimal Strategy (in abs. E) = ', [round(x,3) for x in self.optimal_strategy])

        lin = []
        eff = []

        for i in range(len(self.wells)):

            if abs(self.wells[i].E_lin) < volume_error:
                lin.append(0.)
            else:
                lin.append(round(self.optimal_strategy[i] / self.wells[i].E_lin,3))

            if abs(self.wells[i].E_eff) < volume_error:
                eff.append(0.)
            else:
                eff.append(round(self.optimal_strategy[i] / self.wells[i].E_eff,3))

        print('Optimal Strategy (in E_lin) = ', lin)
        print('Optimal Strategy (in E_eff) = ', eff)
        print('Optimal Strategy (Q) = ', self.Q_wells)
        print('Optimal Strategy (I) = ', self.I_wells)
        print('Optimal Strategy (P) = ', self.P_wells)
        print('I injected = ', sum(self.I_wells))
        print('E injected = ', self.Einj)
        print('P produced = ', sum(self.P_wells))



def RF_opt(x, wells):
    return -sum([float(w.Q_E_interpolate(x[i])) for i, w in enumerate(wells)])


def main():

    ###############################################
    # block for options for OPTIMAL solution search

    write_result_matrix_opt = None      # create txt with wells usability statistics
    display_opt = False                 # statistics of constrained optimization
    multistart_opt = 500                # multistart for constrained optimization problem
    bounds_opt = [0.98]#, 0.95]#, 0.90, 0.85, 0.80] # bounds for E_lin for optimal solution search
    print_log_opt = True                # print current years when profiles are being calculated
    show_graphs_opt = False             # show python figures of profiles after it has been calculated

    # tune some flags
    write_summary_opt = True
    summary_file_opt = os.path.join('./', 'SUMMARY_econ+rf.txt')
    create_schedule_opt = False
    schedule_file_opt = [os.path.join('./', 'GDM_SCHEDULE_'), '.txt']

    ####################################################################################################
    # DATA pre-processing block

    E_annual = [E_mult * p / volume_unit for p in portion]
    E_total = sum(E_annual)

    # get wells from txt files
    wells = GetWellsFromPreparedTxtFiles(data_dir=data_dir, Q_E_str='Q_E_', Q_I_str='Q_I_', P_I_str='P_I_',
                                         well_names=names, well_masks=names_mask, init_option ='files',
                                         E_total=E_total, volume_unit=volume_unit)

    # sort wells: firs well is well with maximum Q(E_lin)
    wells.sort(key=lambda well: max(well.Q_E), reverse=True)

    # define "the most effective point" for each well
    for well in wells:
        well.FindPointOfEffectiviness()

    ####################################################################################################
    # Q MAX block

    report_optimal = Report()
    if write_summary_opt:
        report_optimal.WriteHeader(summary_file_opt)

    for b in bounds_opt:
        solution_naive = SolveConstrainedOptimizationProblem(wells, bounds=[(0, b * w.E_lin) for w in wells],
                                                             multistart=multistart_opt, E_total=E_total,
                                                             write_result_matrix=write_result_matrix_opt,
                                                             display=display_opt)
        solution_naive.PrintInfo()

        # Calculate all sequential strategies and find optimum of economics efficiency.
        # First, construct all sequential strategies
        sequential_strategies = ConstructAllSequentialStrategies(solution_naive.optimal_strategy)

        print('Sequental strategies:')
        print(sequential_strategies)

        # second, find economic efficiency for them

        for ss in sequential_strategies:
            distr = [0.6298, 0, 0, 0.3012, 0, 0, 0, 0, 0, 0, 0, 0.0689, 0, 0]
            assert len(ss) == len(distr)
            ss = [None] * len(distr)
            ss[3], ss[0], ss[11] = 0, 1, 2
            # distr = solution_naive.optimal_strategy
            profiles = CalculateProfiles2(wells, ss, distr, E_annual, first_year,
                                          print_log=print_log_opt, show_graphs=show_graphs_opt,
                                          max_compr_injection_annual=max_compr_injection,
                                          max_compressors_number=Ncompr)

            # calculate economics
            profiles.CalculateEconomics(netback, royalty, income_tax, c_oil, c_gas, r, first_year_econ)

            # append data to report
            if write_summary_opt:
                report_optimal.WriteResults(summary_file_opt, profiles, solution_naive)

            # create schedule file for gdm
            if create_schedule_opt:
                CreateSchedule(schedule_file=schedule_file_opt, profiles=profiles,
                               gdm_mask=gdm_mask, wefac_water=WEFAC, wconinjh_water=WCONINJH)
            break


def SolveConstrainedOptimizationProblem(wells, bounds, multistart, E_total, write_result_matrix, display):

    print('\n#############################################################################')
    print('Constrained Optimization Problem')

    x_opt = []
    f_opt = []

    for m in range(multistart):

        x0 = [random.uniform(0, 1) for w in wells]
        sx = sum(x0)
        for i, xi in enumerate(x0): # TODO (bug was here)
            x0[i] = xi * E_total / sx

        x0 = np.array(x0)
        # print(x0)

        eq_cons = {
            'type': 'eq',
            'fun': lambda x: np.array([sum([x[i] for i in range(len(wells))]) - E_total]),
            'jac': lambda x: np.array([1.0 for i in range(len(wells))])
        }

        res = minimize(RF_opt, x0, method='SLSQP',
                       constraints=[scipy.optimize.LinearConstraint([1 for w in wells], E_total, E_total)],
                       options={'ftol': 1e-8, 'disp': display},
                       bounds=bounds,
                       args=wells)

        if round(-1. * res['fun'], 3) in f_opt:
            pass
        else:
            x_opt.append(res['x'])
            f_opt.append(round(-1. * res['fun'], 3))

    print('End of Constrained Optimization Problem')
    print('#############################################################################\n')

    return OptimalSolution('Constrained Optimization Problem', wells, x_opt[f_opt.index(max(f_opt))], max(f_opt))



def ConstructAllSequentialStrategies(E_list):
    # determine how many wells are involved in injection
    stuff = []
    for e in E_list:
        if abs(e) > volume_error:
            stuff.append(len(stuff))

    # construct all possible sequential strategies for this number of wells
    short_strategies = []
    for i in range(0, len(stuff) + 1):
        for subset in permutations(stuff, i):
            if len(subset) == len(stuff):
                short_strategies.append(subset)

    # add 'None' to non-involved wells
    full_strategies = []
    for st in short_strategies:
        local_counter = 0
        full_strategies.append([])
        for e in E_list:
            if abs(e) < volume_error:
                full_strategies[-1].append(None)
            else:
                full_strategies[-1].append(st[local_counter])
                local_counter += 1

    print('Possible sequential strategies number: ', len(full_strategies))

    return full_strategies




if __name__ == '__main__':
    main()
