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
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.optimize import minimize, rosen, rosen_der
from matplotlib.ticker import FuncFormatter
from itertools import permutations

from chatylka_input_data import *

volume_error = 0.0001

# I
dI = 0.01
I_max = 10

class Well:

    def __init__(self, name, option, data):

        self.name = name

        if option == 'numerical':

            self.S1 = data[0]
            self.A1 = data[1]
            self.I1 = data[2]

            self.A2 = data[3]
            self.I2 = data[4]

            self.data_ = MainDefines(self.S1, self.A1, self.I1, self.A2, self.I2)

            self.I = self.data_[0]
            self.Q_I = self.data_[1]
            self.P_I = self.data_[2]
            self.E = self.data_[3]
            self.Q_E = self.data_[4]

        elif option == 'files':

            self.I = data[0]
            self.Q_I = data[1]
            self.P_I = data[2]
            self.E = data[3]
            self.Q_E = data[4]

            #print(len(self.I), len(self.Q_I), len(self.P_I), len(self.E), len(self.Q_E))

        else:
            print('Error in well definition option\n')
            exit(1)

        self.I_eff = None
        self.Q_eff = None
        self.P_eff = None
        self.E_eff = None
        self.Q_E_eff = None

        self.E_lin = GetElinear(self.E, self.Q_E)
        self.Q_E_interpolate = interp1d(self.E, self.Q_E, 'linear')

        # variables for calculation of the injection dynamics
        self.j = None
        self.dg_dI = None
        self.dH_dE = None

    #
    def Cut(self, E_max, N):

        self.Q_E_interpolate = interp1d(self.E, self.Q_E, 'linear')

        self.E_cut = [ i*E_max/N for i in range(N+1) ]
        self.Q_E_cut = [ self.Q_E_interpolate(e) for e in self.E_cut ]
        self.Q_E_cut_reverse = [self.Q_E_cut[len(self.Q_E_cut) - i - 1] for i in range(len(self.Q_E_cut))]

    def FindPointOfEffectiviness(self):

        eff = 0
        ind_eff = 0

        # go backward till Q < 0 or till i = 1 (because I[0] = 0)
        for i in range(len(self.I) - 1, 0, -1):
            if self.Q_I[i] < 0.:
                break
            else:
                x = self.Q_I[i] / self.I[i]
                if x > eff:
                    eff = x
                    ind_eff = i

        self.I_eff = self.I[ind_eff]
        self.Q_eff = self.Q_I[ind_eff]
        self.P_eff = self.P_I[ind_eff]

        E = self.I_eff - self.P_eff
        E_ = min(self.E, key=lambda y: abs(y - E))

        self.E_eff = E_
        self.Q_E_eff = self.Q_E[(self.E).index(E_)]

    def FindIeffBound(self, bound_constraint):

        I_up = self.I_eff * (1.0 + bound_constraint)
        I_up_closest = min(self.I, key=lambda y: abs(y - I_up))
        I_up_closest_index = (self.I).index(I_up_closest)

        E_up = self.I[I_up_closest_index] - self.P_I[I_up_closest_index]
        E_up_closest = min(self.E, key=lambda y: abs(y - E_up))

        if abs(self.E_eff) < volume_error:
            self.E_bound_mult = 0.
        else:
            self.E_bound_mult = E_up_closest / self.E_eff

    def FindInjectionLimit(self, e, criteria, facility_limit):

        if criteria == 'max':

            E_max = max(e)
            I_max = facility_limit

            dg_rev_dE_max = I_max / E_max
            result_index = 0

            if dg_rev_dE_max < 1.:
                print('Error: too hard condition for max I volume calculation! Exit.')
                exit(1)
            else:

                for i in range(1, len(self.I)):
                    dg_rev_dE_i = (self.I[i] - self.I[i - 1]) / ((self.I[i] - self.P_I[i]) - (self.I[i - 1] - self.P_I[i - 1]))
                    if dg_rev_dE_i > dg_rev_dE_max: # if derivative is bigger than its limit -> stop
                        break
                    else: # go further
                        result_index = i

                self.I_max = self.I[result_index]
                self.E_max = self.I[result_index] - self.P_I[result_index]

        else:
            print('Error: wrong criteria in FindInjectionLimit. Exit')
            exit(1)

    def CalculateDynamicDerivatives(self, E_current): # I = J(E), j = dI/dE

        # find closest point with current E
        E_closest = min(self.E, key=lambda y: abs(y - E_current))
        E_ind = (self.E).index(E_closest)
        Q = self.Q_E[E_ind]
        I_ind = (self.Q_I).index(Q)

        if I_ind == 0:
            self.j = (self.I[I_ind + 1] - self.I[I_ind]) / (
                (self.I[I_ind + 1] - self.P_I[I_ind + 1]) - (self.I[I_ind] - self.P_I[I_ind]))
            self.dg_dI = (self.P_I[I_ind + 1] - self.P_I[I_ind]) / (self.I[I_ind + 1] - self.I[I_ind])
            self.dH_dE = (self.Q_I[I_ind + 1] - self.Q_I[I_ind]) / (
                (self.I[I_ind + 1] - self.P_I[I_ind + 1]) - (self.I[I_ind] - self.P_I[I_ind]))
        else:
            self.j = (self.I[I_ind] - self.I[I_ind - 1]) / (
                (self.I[I_ind] - self.P_I[I_ind]) - (self.I[I_ind - 1] - self.P_I[I_ind - 1]))
            self.dg_dI = (self.P_I[I_ind] - self.P_I[I_ind - 1]) / (self.I[I_ind] - self.I[I_ind - 1])
            self.dH_dE = (self.Q_I[I_ind] - self.Q_I[I_ind - 1]) / (
                (self.I[I_ind] - self.P_I[I_ind]) - (self.I[I_ind - 1] - self.P_I[I_ind - 1]))

        return

class Point:

    def __init__(self, E, wells):
        self.E = E
        for i,e in enumerate(self.E):
            if abs(e) < 0.001:
                self.E[i] = 0.
        self.Q = [ wells[i].Q_E_interpolate(self.E[i]) for i in range(len(wells))]

        self.optimal_value = -1

        self.optimal_pointer = -1

class OptimalSolution:

    def __init__(self, method, wells, strategy, Qmax):

        # check 1
        if len(wells) != len(strategy):
            print('Error: in KnapsackSolution length of wells and strategy lists are not equal!\n')
            exit(1)
        # end check 1

        self.method = method
        self.Qmax = Qmax
        self.optimal_strategy = strategy
        self.wells = wells

        self.Einj = 0.
        for i in range(len(wells)):
            self.Einj += strategy[i]

        self.Q_wells = []
        self.I_wells = []
        self.P_wells = []

        for i,w in enumerate(wells):

            Q = w.Q_E_interpolate(self.optimal_strategy[i])

            Q_closest = min(w.Q_I, key=lambda y: abs(y - Q))
            closest_index = (w.Q_I).index(Q_closest)

            if abs(closest_index - len(w.I)) <= 2:
                print('WARNING: for well ', w.name, ' injected volume = INF., but calculated as finite volume.')
                print('closest index = ', closest_index, ' out of ', len(w.I))

            I = (w.I)[closest_index]
            P = (w.P_I)[closest_index]

            self.Q_wells.append(round(float(Q),3))
            self.I_wells.append(round(I,3))
            self.P_wells.append(round(P,3))

        # check 2
        Qmax_check = sum(self.Q_wells)
        if abs(Qmax_check - self.Qmax) / abs(self.Qmax) > 0.01:
            print('Error: inconsistent data in KnapsackSolution:')
            print('Qmax (external) = ', Qmax, '; Qmax (internal) = ', Qmax_check, '. Exit...')
            exit(1)
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

class FieldDynamics:

    def __init__(self, times, data_tr, data_mer, data_ann, compr_num_float, compr_num, max_compressors_number_at_work, exceed_flags, i_mer_wells):

        self.time = times[0]
        self.dates = times[1]
        self.mer_dates = times[2]
        self.mer_time = times[3]
        self.year_dates = times[4]

        self.q_tr = data_tr[0]
        self.i_tr = data_tr[1]
        self.p_tr = data_tr[2]
        self.e_tr = data_tr[3]

        self.q_mer = data_mer[0]
        self.i_mer = data_mer[1]
        self.p_mer = data_mer[2]
        self.e_mer = data_mer[3]

        self.q_ann = data_ann[0]
        self.i_ann = data_ann[1]
        self.p_ann = data_ann[2]
        self.e_ann = data_ann[3]

        # calculate cumulative volumes
        self.Q_total = sum(self.q_tr)
        self.I_total = sum(self.i_tr)
        self.P_total = sum(self.p_tr)
        self.E_total = sum(self.e_tr)

        self.N_compr = compr_num
        self.N_compr_float = compr_num_float

        self.max_compressors_at_work_total = max_compressors_number_at_work[0]
        self.max_compressors_at_work_uninterruptedly = max_compressors_number_at_work[1]

        self.exceed_1 = exceed_flags[0]
        self.exceed_2 = exceed_flags[1]

        self.i_mer_wells = i_mer_wells

        self.FCF = None
        self.DCF = None
        self.NPV = None
        self.DCP = None

    def CalculateEconomics(self, netback, royalty, income_tax, c_oil, c_gas, r, first_year_econ):

        economics_shift = self.year_dates[0] - first_year_econ

        I = [self.i_ann[i] * volume_unit * Bg for i in range(len(self.year_dates))]
        Q = [self.q_ann[i] * volume_unit * rho / Bo for i in range(len(self.year_dates))]

        self.Revenue = [Q[i] * (netback - royalty - c_oil) for i in range(len(self.year_dates))]
        self.OPEX = [I[i] * c_gas for i in range(len(self.year_dates))]

        self.FCF = [(1.0 - income_tax) * (self.Revenue[i] - self.OPEX[i]) for i in range(len(self.year_dates))]
        self.DCF = [self.FCF[i] / ((1.0 + r) ** (0.5 + i + economics_shift)) for i in range(len(self.year_dates))]
        self.NPV = sum(self.DCF)

        # print(self.Revenue)
        # print(sum(self.Revenue))
        # print(self.OPEX)
        # print(sum(self.OPEX))
        # print(self.FCF)
        # print(sum(self.FCF))
        # print(self.DCF)
        # print(sum(self.DCF))


        # discounted cumulative production (just for interest)
        self.DCP = sum([self.q_ann[i] / ((1 + r) ** (0.5 + i)) for i in range(len(self.year_dates))])

class Report:

    def __init__(self):
        pass

    def WriteHeader(self, result_file):

        f = open(result_file, 'a')
        f.write('NPV\tDCP\tRF_dynamic\tRF_static\tI_dynamic\tI_static\tN_comp1\tN_comp30\tN_comp365\n')
        f.close()

    def WriteResults(self, result_file, profiles, solution_naive):

        print('Writing SUMMARY (NPV, RF, Q, I, P ...) to txt...')
        f = open(result_file, 'a')
        f.write(str(profiles.NPV) + '\t' + str(profiles.DCP) + '\t' + str(profiles.Q_total) + '\t' +
                str(solution_naive.Qmax) + '\t' + str(profiles.I_total) + '\t' + str(sum(solution_naive.I_wells))
                + '\t' + str(profiles.N_compr_float[0]) + '\t' + str(profiles.N_compr_float[1]) + '\t' +
                str(profiles.N_compr_float[2]) + '\t' + str(profiles.max_compressors_at_work_total) +
                '\t' + str(profiles.max_compressors_at_work_uninterruptedly) +
                '\t' + str(profiles.exceed_1) + '\t' + str(profiles.exceed_2) + '\n')
        f.close()

def RF_opt(x, wells):
    Q = 0.
    for i,w in enumerate(wells):
        Q += w.Q_E_interpolate(x[i])
    #print(Q)
    return -1.*Q

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
    summary_file_opt = os.path.join(data_dir, 'SUMMARY_econ+rf.txt')
    create_schedule_opt = True
    schedule_file_opt = [os.path.join(data_dir, 'GDM_SCHEDULE_'), '.txt']

    ###############################################
    # block for options for BASE solution search

    base_case = True           # calculate base case?

    write_result_matrix_base = None  # create txt with wells usability statistics
    display_base = False             # statistics of constrained optimization
    multistart_base = 500            # multistart for constrained optimization problem
    bounds_base = [0.30]#[0.05, 0.1, 0.15, 0.20, 0.25, 0.30]  # bounds for I_eff for optimal solution search
    print_log_base = True            # print current years when profiles are being calculated
    show_graphs_base = True         # show python figures of profiles after it has been calculated

    # tune some flags
    write_summary_base = True
    summary_file_base = os.path.join(data_dir, 'SUMMARY_BASE_econ+rf.txt')
    create_schedule_base = True
    schedule_file_base = [os.path.join(data_dir, 'GDM_SCHEDULE_BASE_'), '.txt']

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

            # get profiles
            #
            profiles = CalculateProfiles2(wells, ss, solution_naive.optimal_strategy, E_annual, first_year,
                                             print_log=print_log_opt, show_graphs=show_graphs_opt,
                                          max_compr_injection_annual=max_compr_injection,
                                          max_compressors_number=Ncompr)

            # check cumulative volumes
            CheckCumulativeVolumes(profiles, solution_naive)

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
    ####################################################################################################
    # Base approach for Q MAX

    if base_case:

        report_base = Report()
        if write_summary_base:
            report_base.WriteHeader(summary_file_base)

        # # solve Knapsack problem. The disadvantage of this method - E_total >= E_injected
        # solution1 = KnapsackSolutionBruteForce(wells=wells, capacity=E_total, make_add=False)
        # solution1.PrintInfo()
        #
        # # another solution - Knapsack with additional injection
        # solution2 = KnapsackSolutionBruteForce(wells = wells, capacity = E_total, make_add = True)
        # solution2.PrintInfo()

        # solve Knapsack problem with accuracy +- delta. That is equivalent to ordirar constrained optimization
        # bound_constraint = 0.3 # bounds_I = [0, 1 + bound_constraint * I_eff]
        for b in bounds_base:

            for well in wells:
                well.FindIeffBound(b)
            solution_pseudo_knapsack = SolveConstrainedOptimizationProblem(wells, bounds=[(0.0, w.E_bound_mult * w.E_eff) for w in wells],
                                                                 multistart=multistart_base, E_total=E_total,
                                                                 write_result_matrix=write_result_matrix_base,
                                                                           display=display_base)
            solution_pseudo_knapsack.PrintInfo()

            sequential_strategies = ConstructAllSequentialStrategies(solution_pseudo_knapsack.optimal_strategy)

            print('Sequental strategies:')
            print(sequential_strategies)

            for ss in sequential_strategies:

                # get profiles
                profiles = CalculateProfiles2(wells, [2, None, None, 3, 0, None, 1, None, None, None, None, None, None, None], solution_pseudo_knapsack.optimal_strategy, E_annual, first_year,
                                             print_log=print_log_base, show_graphs=show_graphs_base,
                                              max_compr_injection_annual=max_compr_injection,
                                              max_compressors_number=Ncompr)

                # check cumulative volumes
                CheckCumulativeVolumes(profiles, solution_pseudo_knapsack)

                # calculate economics
                profiles.CalculateEconomics(netback, royalty, income_tax, c_oil, c_gas, r, first_year_econ)

                # append data to report
                if write_summary_base:
                    report_base.WriteResults(summary_file_base, profiles, solution_pseudo_knapsack)

                # create schedule file for gdm
                if create_schedule_base:
                    CreateSchedule(schedule_file=schedule_file_base, profiles=profiles,
                                   gdm_mask=gdm_mask, wefac_water=WEFAC, wconinjh_water=WCONINJH)
                break

    # Optimization: brutforce or dynprog
    #opt_str = BrutForce(n_strategies = 20, annual_injection = E_annual, discount = r,
    #          start_point=np.zeros(2), end_point=np.array([E_opt1, E_opt2]), wells = wells)

    # DCP optimization
    #DCPOptimization(n_points=17, annual_injection=E_annual, discount=r,
    #          start_point=np.zeros(2), wells=wells)

    #DynProg(E_annual, YearsNum, [E_opt1, E_opt2], wells, r, 40)

################################################################################
################################################################################
def CreateSchedule(schedule_file, profiles, gdm_mask, wefac_water, wconinjh_water):

    year_grp_ctrl_str = 'WUCTRL'
    well_ctrl_str = 'WUC'
    result = ''
    # flag = 'water': water is injected and no gas was injected previously
    # flag = 'gas': gas is injected
    # flag = 'water+': water is injected after gas injection
    new_well_for_inj = ['water' for w in range(len(profiles.i_mer_wells[0]))]

    # time loop
    for t, date in enumerate(profiles.mer_dates):

        # 0) define the current year and current month
        current_year = (eval(date.split('.')[1]) - eval((profiles.mer_dates[0]).split('.')[1])) + 1
        current_month = eval(date.split('.')[0])

        # 1) UDQ block
        udq_str = 'UDQ\n'

        total_rate = sum([profiles.i_mer_wells[w + 1][t] for w in range(len(profiles.i_mer_wells[0]))])

        for w in range(len(profiles.i_mer_wells[0])):

            well_portion = profiles.i_mer_wells[w + 1][t] / total_rate

            # write UDQ variable only if rate != 0 and well injects gas first time
            if round(well_portion,2) > 0.:
                udq_str += 'DEFINE\t' + well_ctrl_str + str(profiles.i_mer_wells[0][w]) + '\t' \
                           + str(round(well_portion,2)) + ' * ' + year_grp_ctrl_str + str(current_year) \
                           + '\t/\n'
                new_well_for_inj[w] = 'gas'
            elif round(well_portion,2) == 0. and new_well_for_inj[w] == 'gas':
                # special case: injector starts to inject water again
                new_well_for_inj[w] = 'water+'
        udq_str += '/\n'

        result += udq_str + '\n'

        # 2) DATES block
        # form 01.2021 make 1 JAN 2021
        date_str = ConvertDate(date)
        result += date_str + '\n'

        # 3) GCONPROD and GCONINJE blocks
        if current_month == 1:
            gcon_str = 'GCONPROD\nFIELD1 GRAT 2* GUMAXGP' + str(current_year) + ' 1* TARG NO /\n/\n'
            gcon_str += 'GCONINJE\nINJWELLS GAS RATE GUGIR' + str(current_year) + ' /\n/\n'
            result += gcon_str + '\n'

        # 4) WEFAC block for gas injectors
        wefac_str = 'WEFAC\n'
        for w in range(len(profiles.i_mer_wells[0])):
            if new_well_for_inj[w] == 'gas':
                gdm_well_name = gdm_mask[profiles.i_mer_wells[0][w]]
                wefac_str += "'" + str(gdm_well_name) + "'" + '\t1.0' + '\t/\n'
        wefac_str += '/\n'
        result += wefac_str + '\n'

        # 5) WCONINJE block for gas injectors
        wconinje_str = 'WCONINJE\n'
        for w in range(len(profiles.i_mer_wells[0])):
            if new_well_for_inj[w] == 'gas':
                gdm_well_name = gdm_mask[profiles.i_mer_wells[0][w]]
                wconinje_str += "'" + str(gdm_well_name) + "'" + '\tGAS\tOPEN\tRATE\t' \
                                + well_ctrl_str + str(profiles.i_mer_wells[0][w]) + '\t/\n'
        wconinje_str += '/\n'
        result += wconinje_str + '\n'

        # 6) WEFAC block for wells which started to inject water AGAIN
        wefac_str = 'WEFAC\n'
        for w in range(len(profiles.i_mer_wells[0])):
            if new_well_for_inj[w] == 'water+':
                gdm_well_name = gdm_mask[profiles.i_mer_wells[0][w]]
                wefac_str += "'" + str(gdm_well_name) + "'" + '\t' + str(wefac_water[profiles.i_mer_wells[0][w]]) \
                             + '\t/\n'
        wefac_str += '/\n'
        result += wefac_str + '\n'

        # 7) WCONINJH block for wells which started to inject water AGAIN
        wconinjh_str = 'WCONINJE\n'
        for w in range(len(profiles.i_mer_wells[0])):
            if new_well_for_inj[w] == 'water+':
                gdm_well_name = gdm_mask[profiles.i_mer_wells[0][w]]
                wconinjh_str += "'" + str(gdm_well_name) + "'" + '\t' + 'WATER OPEN	RATE\t' \
                               + str(wconinjh_water[profiles.i_mer_wells[0][w]]) + '\t/\n'
        wconinjh_str += '/\n'
        result += wconinjh_str + '\n'

        # special correction: don't repeat water injection keywords next step
        for w in range(len(profiles.i_mer_wells[0])):
            if new_well_for_inj[w] == 'water+':
                new_well_for_inj[w] = 'water'

    # write the last date
    last_date = IncrementMonthData(profiles.mer_dates[-1])
    date_str = ConvertDate(last_date)
    result += date_str + '\n'

    # write to file
    wells_str = ''
    for w in range(len(profiles.i_mer_wells[0])):
        wells_str += 'w' + str(profiles.i_mer_wells[0][w]) + '_'
    my_file = schedule_file[0] + wells_str + str(round(profiles.Q_total, 2)) + '_' \
              + str(round(profiles.NPV, 2)) + '_' + schedule_file[-1]
    f = open(my_file, 'w')
    f.write(str(result))
    f.close()

    return

def IncrementMonthData(date):

    m, y = date.split('.')
    m = eval(m)
    y = eval(y)

    if m == 12:
        y += 1
        m = 1
    else:
        m += 1

    next_date = str(m) + '.' + str(y)

    return next_date

def ConvertDate(date):

    months = {
        1: 'JAN',
        2: 'FEB',
        3: 'MAR',
        4: 'APR',
        5: 'MAY',
        6: 'JUN',
        7: 'JUL',
        8: 'AUG',
        9: 'SEP',
        10: 'OCT',
        11: 'NOV',
        12: 'DEC'
    }

    month_str, year_str = date.split('.')
    month_str = months[eval(month_str)]

    return('DATES\n' + '1\t' + month_str + '\t' + year_str + '\t' + '/\n/\n')

def GetWellsFromPreparedTxtFiles(data_dir, Q_E_str, Q_I_str, P_I_str,
                                         well_names, well_masks, init_option, E_total, volume_unit):

    wells = []

    for i in range(len(well_names)):
        f = open(os.path.join(data_dir, Q_E_str + str(well_names[i]) + '.txt'), 'r')
        lines = f.readlines()
        f.close()
        E = []
        Q = []
        for line in lines:
            EQ = line.split('\t')
            E.append(eval(EQ[0]) / volume_unit)
            Q.append(eval(EQ[1]) / volume_unit)

        f = open(os.path.join(data_dir, Q_I_str + str(well_names[i]) + '.txt'), 'r')
        lines = f.readlines()
        f.close()
        I = []
        Q_I = []
        for line in lines:
            IQ = line.split('\t')
            I.append(eval(IQ[0]) / volume_unit)
            Q_I.append(eval(IQ[1]) / volume_unit)

        f = open(os.path.join(data_dir, Q_I_str + str(well_names[i]) + '.txt'), 'r')
        lines = f.readlines()
        f.close()
        I_ = []
        P_I = []
        for line in lines:
            IP = line.split('\t')
            I_.append(eval(IP[0]) / volume_unit)
            P_I.append(eval(IP[1]) / volume_unit)

        if I != I_:
            print('Error in file loading for well ' + str(well_names[i]) + '\n')
            exit(1)

        wells.append(Well(well_masks[str(well_names[i])], init_option, [I, Q_I, P_I, E, Q]))

    # check the loaded data
    E_lin_summ = 0.
    for i, w in enumerate(wells):
        E_lin_summ += w.E_lin
        # print("E_lin ", str(w.name), " = ", w.E_lin)
    if E_total >= E_lin_summ:
        print('Error: trivial solution: E = E_lin')
        print('E_total = ', E_total, '; but the summ of E_lin = ', E_lin_summ)
        print('Exit...')
        exit(1)

    return wells

def SolveConstrainedOptimizationProblem(wells, bounds, multistart, E_total, write_result_matrix, display):

    print('\n#############################################################################')
    print('Constrained Optimization Problem')

    x_opt = []
    f_opt = []

    for m in range(multistart):

        x0 = [random.uniform(0, 1) for w in wells]
        sx = sum(x0)
        for xi in x0:
            xi *= E_total / sx

        x0 = np.array(x0)
        # print(x0)

        eq_cons = {
            'type': 'eq',
            'fun': lambda x: np.array([sum([x[i] for i in range(len(wells))]) - E_total]),
            'jac': lambda x: np.array([1.0 for i in range(len(wells))])
        }

        res = minimize(RF_opt, x0, method='SLSQP',
                       constraints=[eq_cons],
                       options={'ftol': 1e-8, 'disp': display},
                       bounds=bounds,
                       args=wells)

        # print(res)
        if round(-1. * res['fun'], 3) in f_opt:
            pass
        else:
            #print(m + 1, -1. * res['fun'])
            x_opt.append(res['x'])
            f_opt.append(round(-1. * res['fun'], 3))

    Q_opt = []
    for i, w in enumerate(wells):
        Q_opt.append(float(w.Q_E_interpolate(x_opt[f_opt.index(max(f_opt))][i])))
        if Q_opt[-1] < 0.001:
            Q_opt[-1] = 0.0

    if write_result_matrix != None:

        print('Write results to txt file...')

        data_dir = write_result_matrix[0]
        filename = write_result_matrix[1]

        ff = copy.copy(f_opt)
        ff.sort()
        xx = []
        for f in ff:
            xx.append(x_opt[f_opt.index(f)])

        f = open(data_dir + '\\' + filename, 'w')
        for i,fi in enumerate(ff):
            f.write(str(fi))
            for j,x in enumerate(xx[i]):
                if abs(x) > volume_error:
                    f.write('\t' + str(j + 1))
                else:
                    f.write('\t' + '0')
            f.write('\n')
        f.close()

    print('End of Constrained Optimization Problem')
    print('#############################################################################\n')

    return OptimalSolution('Constrained Optimization Problem', wells, x_opt[f_opt.index(max(f_opt))], max(f_opt))

def KnapsackSolutionBruteForce(wells, capacity, make_add):

    print('\n#############################################################################')
    print('Knapsack Optimization Problem')

    # 0. print info
    # print(capacity)
    # for well in wells:
    #     print(well.E_eff, well.E_lin, well.E_eff / well.E_lin, ' ---- ', well.Q_eff, well.Q_E[-1], well.Q_eff / well.Q_E[-1])

    # 1. create all knapsack strategies: N = 2^len(wells)
    knapsack_strategies = [[1], [0]]
    for i in range(len(wells) - 1): #(3)#
        knapsack_strategies = AppendStrategyToTree(knapsack_strategies)
    print('All knapsack strategies number = ', len(knapsack_strategies))
    #print(knapsack_strategies)

    # 2. filter strategies according to capacity and calculate Q
    knapsack_strategies_UPD = []
    Q = []

    min_E_eff = min(well.E_eff for well in wells if well.E_eff > 0)

    for strategy in knapsack_strategies:

        weight = 0.
        for i,well in enumerate(wells):
            weight += well.E_eff * strategy[i]

        if weight <= capacity:
            knapsack_strategies_UPD.append(strategy)
            Qstr = 0.
            for j, well in enumerate(wells):
                Qstr += well.Q_eff * strategy[j]

            # if weight < capacity, we can inject to at least one well
            dw = capacity - weight
            dQ = []
            if dw <= min_E_eff and dw > 0. and make_add == True:
                for j, well in enumerate(wells):
                    if strategy[j] == 0:
                        # try to inject to a new well
                        dQ.append(well.Q_E_interpolate(dw))
                    else:
                        # try to inject additional volume to involved well
                        dQ.append(well.Q_E_interpolate(well.E_eff + dw) - well.Q_eff)
            else:
                dQ.append(0.)
            add_inj = dQ.index(max(dQ))
            if max(dQ) > 0:
                strategy[add_inj] += dw / wells[add_inj].E_eff

            Q.append(Qstr + max(max(dQ),0))
    print('Knapsack strategies number with propper weight = ', len(knapsack_strategies_UPD))

    # convert E_eff to abs.vol. for optimal strategy
    opt_str = copy.copy(knapsack_strategies_UPD[Q.index(max(Q))])

    print('End of Knapsack Optimization Problem')
    print('#############################################################################\n')

    return OptimalSolution('Knapsack', wells, [x * wells[j].E_eff for j,x in enumerate(opt_str)], max(Q))

def AppendStrategyToTree(knapsack_strategies):

    children = []

    for i,strategy in enumerate(knapsack_strategies):
        c1 = copy.copy(strategy)
        c2 = copy.copy(strategy)
        c1.append(1)
        c2.append(0)
        children.append(c1)
        children.append(c2)

    return children

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
        full_strategies.append([])
        local_counter = 0
        for e in E_list:
            if abs(e) < volume_error:
                full_strategies[-1].append(None)
            else:
                full_strategies[-1].append(st[local_counter])
                local_counter += 1

    print('Possible sequential strategies number: ', len(full_strategies))

    return full_strategies

def CalculateProfiles2(wells_all, sequence, E_distribution_all, e_annual, first_year, print_log, show_graphs, max_compr_injection_annual, max_compressors_number):

    bug_mode = False

    print('Calculating dynamics...')
    # now we consider only sequential modes
    # sequence = [None, None, 0, 1, None, None, 2, None] - it means that only 3 wells are involved and numbers mean the sequence

    if print_log:
        print('sequence = ', sequence)

    # resulted profiles
    q_field = []
    i_field = []
    p_field = []
    e_field = []

    if len(wells_all) != len(sequence):
        print('Error in lengths in CalculateProfiles. Exit...')
        exit(1)

    # daily injection limit
    max_i_daily = max_compr_injection_annual * max_compressors_number / volume_unit / 365.

    # how many wells we should inject
    wells_num = 0
    for i in range(len(wells_all)):
        if type(sequence[i]) == int:
            wells_num += 1
    # local variables determination
    wells = [0. for i in range(wells_num)] # wells for injection
    E_distribution = [0. for i in range(wells_num)] # final volumes of external supplies to each well
    wells_injected_external_volume = [0. for i in range(wells_num)] # current volumes of external supplies to each well
    mers_for_gdm = [[] for i in range(wells_num)]
    for i, w in enumerate(wells_all):
        if type(sequence[i]) == int:
            wells[sequence[i]] = copy.copy(w)
            E_distribution[sequence[i]] = copy.copy(E_distribution_all[i])
            wells_injected_external_volume[sequence[i]] = 0.
            # calculate current derivatives
            (wells[sequence[i]]).CalculateDynamicDerivatives(0.)

    # current well for injection
    order = 0

    # flag means that injection has exceeded the limit for 1,2,...N-1 wells
    exceed_1 = False

    # flag means that injection has exceeded the limit for N-th well
    exceed_2 = False

    # year-loop
    for t in range(len(e_annual)):

        # define the current year
        year = first_year + t
        days_in_current_year = (datetime.date(year + 1, 1, 1) - datetime.date(year, 1, 1)).days
        if print_log:
            print('Year ', year, ' days ', days_in_current_year)

        # external supply  rate during the year
        e_daily = e_annual[t] / days_in_current_year

        # day-loop during the year
        for day in range(1, days_in_current_year + 1):

            if bug_mode:
                print(day, ' main well ', (wells[order]).name)

            # We know the current E for the well - wells_injected_external_volume[well_index].
            # Let's find the corresponding I, Q and P values and then find i, q and p
            i_daily = (wells[order]).j * e_daily

            # current injection i_daily can be more than central facility constraint
            if i_daily < max_i_daily:

                if bug_mode:
                    print('no assistant')

                # it's ok: calculate current p and q
                p_daily = (wells[order]).dg_dI * i_daily
                q_daily = (wells[order]).dH_dE * e_daily

                # take into account the injected supplies
                wells_injected_external_volume[order] += e_daily * 1.0

                # recalculate the dynamics properties of the well
                (wells[order]).CalculateDynamicDerivatives(wells_injected_external_volume[order])

                for w in range(wells_num):
                    if w == order:
                        mers_for_gdm[w].append(i_daily)
                    else:
                        mers_for_gdm[w].append(0.)

            else:
                # we should inject some volume to the next well
                if order < len(wells) - 1:

                    # flag: have we found another well for injection?
                    assistant_well = False

                    # find next well with j < j of the current well
                    for ss in range(1, len(wells) - order):

                        if bug_mode:
                            print('Try ', (wells[order + ss]).name, ' with derivative = ', (wells[order + ss]).j, ' Current derivative = ', (wells[order]).j)

                        # check that portions of external supplies will > 0 and < e_daily
                        e_1 = (max_i_daily - (wells[order + ss]).j * e_daily) / (
                        (wells[order]).j - (wells[order + ss]).j)

                        if ((wells[order + ss]).j < (wells[order]).j) and (e_1 > 0.) and (e_1 < e_daily):

                            if bug_mode:
                                print('assistant is ', wells[order + ss].name)

                            assistant_well = True

                            # portion of external supplies to the next well
                            e_2 = e_daily - e_1

                            if bug_mode:
                                print('e_daily ', e_daily, ' = ', e_2, ' + ', e_1)

                            # calculate rates
                            i_daily = (wells[order]).j * e_1 + (wells[order + ss]).j * e_2
                            p_daily = (wells[order]).dg_dI * (wells[order]).j * e_1 + (wells[order + ss]).dg_dI * (
                                wells[order + ss]).j * e_2
                            q_daily = (wells[order]).dH_dE * e_1 + (wells[order + ss]).dH_dE * e_2

                            # take into account the injected supplies
                            wells_injected_external_volume[order] += e_1 * 1.0
                            wells_injected_external_volume[order + ss] += e_2 * 1.0

                            # recalculate the dynamics properties of the well
                            (wells[order]).CalculateDynamicDerivatives(wells_injected_external_volume[order])
                            (wells[order + ss]).CalculateDynamicDerivatives(wells_injected_external_volume[order + ss])

                            for w in range(wells_num):
                                if w == order:
                                    mers_for_gdm[w].append((wells[order]).j * e_1)
                                elif w == order + ss:
                                    mers_for_gdm[w].append((wells[order + ss]).j * e_2)
                                else:
                                    mers_for_gdm[w].append(0.)

                            break

                    if not assistant_well: # there were no additional well for injection

                        exceed_1 = True

                        if print_log:
                            print(
                                "WARNING: there were no additional well for injection; it exceeds the limit; all is injected to the current well")
                            if bug_mode:
                                input()

                        # it's NOT ok, but we can't do anything: calculate current p and q
                        p_daily = (wells[order]).dg_dI * i_daily
                        q_daily = (wells[order]).dH_dE * e_daily

                        # take into account the injected supplies
                        wells_injected_external_volume[order] += e_daily * 1.0

                        # recalculate the dynamics properties of the well
                        (wells[order]).CalculateDynamicDerivatives(wells_injected_external_volume[order])

                        for w in range(wells_num):
                            if w == order:
                                mers_for_gdm[w].append(i_daily)
                            else:
                                mers_for_gdm[w].append(0.)

                else:

                    exceed_2 = True

                    if print_log:
                        print("WARNING: injection value to the last well is exceeds the limit! But we can't do anything")
                        if bug_mode:
                            input()

                    # it's NOT ok, but we can't do anything: calculate current p and q
                    p_daily = (wells[order]).dg_dI * i_daily
                    q_daily = (wells[order]).dH_dE * e_daily

                    # take into account the injected supplies
                    wells_injected_external_volume[order] += e_daily * 1.0

                    # recalculate the dynamics properties of the well
                    (wells[order]).CalculateDynamicDerivatives(wells_injected_external_volume[order])

                    for w in range(wells_num):
                        if w == order:
                            mers_for_gdm[w].append(i_daily)
                        else:
                            mers_for_gdm[w].append(0.)

            # if injected volume is bigger than E_distribution -> inject to the next well
            if wells_injected_external_volume[order] >= E_distribution[order]:# or abs(E_distribution[order] - wells_injected_external_volume[order]) < volume_error:

                # go to the next well
                order += 1

                if bug_mode:
                    print(wells_injected_external_volume, E_distribution)
                    print('GO TO NEXT WELL\n')
                    input()

                # if all wells were involved, but there are some days in the current year
                if order == len(wells) and day != days_in_current_year:

                    if print_log:
                        print("Warning: there are ", (days_in_current_year - day), " days rest, but all wells were injected! Make all rates = 0.")
                    # exit(1)

                    order -= 1 # now order value doesn't matter, but order must be in list
                    i_daily = 0.
                    p_daily = 0.
                    q_daily = 0.
                    e_daily = 0. # it's important!

            # finally, fill the daily profiles
            q_field.append(q_daily * 1.0)
            p_field.append(p_daily * 1.0)
            i_field.append(i_daily * 1.0)
            e_field.append(e_daily * 1.0)

        if bug_mode:
            input()

    # make monthly reports and annual volumes
    print('Creating profiles...')

    time = [] # 1,2,3,...7305
    dates = [] # string type: '01.01.2021', '02.01.2021', ...
    mer_dates = [] # string type: '01.2021', '02.2021', ...
    mer_time = [] # 1,2,... 240
    year_dates = [] # 2021, 2022, ...

    # FIELD
    q_mer = []
    q_ann = [0.]
    i_mer = []
    i_ann = [0.]
    p_mer = []
    p_ann = [0.]
    e_mer = []
    e_ann = [0.]

    # WELLS
    i_mer_wells = [[] for i in range(wells_num)]

    # local summators
    month_duration = 0.
    sum_q = 0.
    sum_i = 0.
    sum_p = 0.
    sum_e = 0.
    sum_i_wells = [0. for i in range(wells_num)]

    # cut MER and annual profiles from cumulative profiles
    current_date = datetime.date(first_year, 1, 1)  # datetime.timedelta(days=1)

    year_dates.append(current_date.year)

    for t in range(len(e_field)):

        time.append(t + 1)
        dates.append(str(current_date.day) + '.' + str(current_date.month) + '.' + str(current_date.year))

        sum_q += q_field[t]
        sum_i += i_field[t]
        sum_p += p_field[t]
        sum_e += e_field[t]

        for w in range(wells_num):
            sum_i_wells[w] += mers_for_gdm[w][t]

        month_duration += 1

        # if this day is a final day of the month
        if (current_date + datetime.timedelta(days=1)).month != current_date.month:

            # fix current month date
            mer_dates.append(str(current_date.month) + '.' + str(current_date.year))
            mer_time.append(len(mer_time) + 1)

            # fix average rates
            q_mer.append(sum_q / month_duration)
            i_mer.append(sum_i / month_duration)
            p_mer.append(sum_p / month_duration)
            e_mer.append(sum_e / month_duration)

            for w in range(wells_num):
                i_mer_wells[w].append(sum_i_wells[w] / month_duration)

            # add cumulative volumes to the current annual production and injection
            q_ann[-1] += sum_q
            i_ann[-1] += sum_i
            p_ann[-1] += sum_p
            e_ann[-1] += sum_e

            # cut the annual cumulative volumes
            if (current_date + datetime.timedelta(days=1)).year != current_date.year and t != (len(e_field) - 1):
                q_ann.append(0.)
                i_ann.append(0.)
                p_ann.append(0.)
                e_ann.append(0.)
                year_dates.append(current_date.year + 1)

            # reset summators for the next month
            sum_q = 0.
            sum_i = 0.
            sum_p = 0.
            sum_e = 0.
            sum_i_wells = [0. for i in range(wells_num)]
            month_duration = 0.

        # next step
        current_date += datetime.timedelta(days=1)

    # check FIELD and WELLS rates
    if print_log:
        print('Compare field rates and sum of field rates...')
    max_delta = 0.0
    for month in range(len(i_mer)):
        local = sum([ i_mer_wells[w][month] for w in range(wells_num) ])
        if ( abs(local - i_mer[month]) / i_mer[month] * 100. ) > max_delta:
            max_delta = abs(local - i_mer[month]) / i_mer[month] * 100.
            if max_delta > 5.0:
                print('Error: the sum of well rates is not equal to field rate with precision 5%! Exit...')
                exit(1)
    if print_log:
        print('Max internal delta for rates in % = ', max_delta)

    # convert to abs. units: SM3/day
    for w in range(len(i_mer_wells)):
        for t in range(len(i_mer_wells[w])):
            i_mer_wells[w][t] *= volume_unit * Bg

    # add well names
    i_mer_wells.insert(0, [])
    # correct names according to mask
    for w in wells:
        for real_name, mask_name in names_mask.items():
            if mask_name == w.name:
                i_mer_wells[0].append(real_name)

    # calculate how many compressors do we need using different stepsizes: day, month and year
    one_compr_max = [max_compr_injection_annual / volume_unit / 365., max_compr_injection_annual / volume_unit / 365., max_compr_injection_annual / volume_unit ] # in relative units
    compr_num_float = [max(i_field) / one_compr_max[0], max(i_mer) / one_compr_max[1], max(i_ann) / one_compr_max[2]]
    compr_num = [int(compr_num_float[0]), int(compr_num_float[1]), int(compr_num_float[2])]
    for cn in range(len(compr_num)):
        if compr_num_float[cn] != compr_num[cn]:
            compr_num[cn] += 1

    print('Compressors number = ', compr_num_float, ' ', compr_num)

    # calculate how many days TOTALLY and UNINTERRUPTLY do we need the maximum number of working compressors
    if max_compressors_number == 1:
        max_compressors_number_at_work = [len(i_field) / 365., len(i_field) / 365.]
    else:
        # the power of (N-1) compressors
        pre_limit = one_compr_max[0] * (max_compressors_number - 1)
        # periods of injection when we needed the max number of compressors
        periods = [0]
        for k, inj in enumerate(i_field):
            if inj > pre_limit:
                periods[-1] += 1
                # stop the current period if next Inj rate is less than the limit; so we start the next period
                if k != (len(i_field) - 1) and i_field[k + 1] <= pre_limit:
                    periods.append(0)
        max_compressors_number_at_work = [sum(periods) / 365., max(periods) / 365.]
    print('Max compressors at work (in years) = ', max_compressors_number_at_work[0], ' total and', max_compressors_number_at_work[1], ' max constant period')

    if show_graphs:

        plt.plot(time, q_field, color='r')
        plt.plot(time, i_field, color='b')
        plt.plot(time, p_field, color='cyan')
        plt.plot(time, e_field, color='g')
        plt.plot([time[0], time[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[0] + 1):
            plt.plot([time[0], time[-1]], [k * one_compr_max[0], k * one_compr_max[0]], color='maroon', linestyle='-')
        plt.show()

        plt.plot(mer_time, q_mer, color='r', marker='.')
        plt.plot(mer_time, i_mer, color='b', marker='.')
        plt.plot(mer_time, p_mer, color='cyan', marker='.')
        plt.plot(mer_time, e_mer, color='g', marker='.')
        plt.plot([mer_time[0], mer_time[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[1] + 1):
            plt.plot([mer_time[0], mer_time[-1]], [k * one_compr_max[1], k * one_compr_max[1]], color='maroon', linestyle='-')
        plt.show()

        plt.plot(year_dates, q_ann, color='r', marker='.')
        plt.plot(year_dates, i_ann, color='b', marker='.')
        plt.plot(year_dates, p_ann, color='cyan', marker='.')
        plt.plot(year_dates, e_ann, color='g', marker='.')
        plt.plot([year_dates[0], year_dates[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[2] + 1):
            plt.plot([year_dates[0], year_dates[-1]], [k * one_compr_max[2], k * one_compr_max[2]], color='maroon', linestyle='-')
        plt.show()

    # construct the solution
    result = FieldDynamics([time, dates, mer_dates, mer_time, year_dates], [q_field, i_field, p_field, e_field],
                           [q_mer, i_mer, p_mer, e_mer], [q_ann, i_ann, p_ann, e_ann], compr_num_float, compr_num, max_compressors_number_at_work, [exceed_1, exceed_2], i_mer_wells)

    if print_log:
        print(len(i_field))

    return result

def FindRates(well, E_current, e_daily):
    # find closest point with current E
    E_closest = min(well.E, key=lambda y: abs(y - E_current))
    E_ind = (well.E).index(E_closest)
    Q = well.Q_E[E_ind]
    I_ind = (well.Q_I).index(Q)

    if I_ind == 0:
        dg_rev_dE = (well.I[I_ind + 1] - well.I[I_ind]) / (
            (well.I[I_ind + 1] - well.P_I[I_ind + 1]) - (well.I[I_ind] - well.P_I[I_ind]))
        dg_dI = (well.P_I[I_ind + 1] - well.P_I[I_ind]) / (well.I[I_ind + 1] - well.I[I_ind])
        dH_dE = (well.Q_I[I_ind + 1] - well.Q_I[I_ind]) / (
            (well.I[I_ind + 1] - well.P_I[I_ind + 1]) - (well.I[I_ind] - well.P_I[I_ind]))
    else:
        dg_rev_dE = (well.I[I_ind] - well.I[I_ind - 1]) / (
            (well.I[I_ind] - well.P_I[I_ind]) - (well.I[I_ind - 1] - well.P_I[I_ind - 1]))
        dg_dI = (well.P_I[I_ind] - well.P_I[I_ind - 1]) / (well.I[I_ind] - well.I[I_ind - 1])
        dH_dE = (well.Q_I[I_ind] - well.Q_I[I_ind - 1]) / (
            (well.I[I_ind] - well.P_I[I_ind]) - (well.I[I_ind - 1] - well.P_I[I_ind - 1]))

    i_daily = dg_rev_dE * e_daily
    p_daily = dg_dI * i_daily
    q_daily = dH_dE * e_daily

    return i_daily, p_daily, q_daily

def CalculateProfiles(wells, sequence, E_distribution, e_annual, first_year, print_log, show_graphs, inj_annual_constraint):

    print('Calculating dynamics...')
    # now we consider only sequential modes
    # sequence = [None, None, 0, 1, None, None, 2, None] - it means that only 3 wells are involved and numbers mean the sequence

    # resulted profiles
    q_field = []
    i_field = []
    p_field = []
    e_field = []

    # local variable to determine current injected volume to each well
    wells_injected_external_volume = [0. for w in wells]

    if len(wells) != len(sequence):
        print('Error in lengths in CalculateProfiles. Exit...')
        exit(1)

    for t in range(len(e_annual)):

        # define the current year
        year = first_year + t
        days_in_current_year = (datetime.date(year + 1, 1, 1) - datetime.date(year, 1, 1)).days
        if print_log:
            print('Year ', year, ' days ', days_in_current_year)

        q = []
        i = []
        p = []
        e = []

        # # number of periods for integration
        # num_periods = int(days_in_current_year / precision_internal)
        # if days_in_current_year % precision_internal != 0:
        #     num_periods += 1
        #
        # # number of days in each period of integration
        # days_in_periods = [precision_internal for k in range(num_periods - 1)]
        # days_in_periods.append(days_in_current_year - precision_internal * (num_periods - 1))

        # how much gas of the current year was NOT injected
        rest_gas = e_annual[t]

        # external supply  rate during the year
        e_daily = e_annual[t] / days_in_current_year

        # global counter of time in this year; type = float
        current_t = 0.0

        # In year "t" we have injection volume e[t].
        # Let's inject this volume to all wells for which sequence[i] != None.
        # Use wells_injected_external_volume to determine how much volume is possible to inject to each well
        for order in range(len(wells)):

            # search well which is involved in injection
            well_index = sequence.index(order)
            well = wells[well_index]

            # if injected volume is less than E_distribution -> inject to this well
            if wells_injected_external_volume[well_index] < E_distribution[well_index] and abs(E_distribution[well_index] - wells_injected_external_volume[well_index]) > volume_error:

                # Inject gas to this well.
                # How much gas will be injected?
                if (wells_injected_external_volume[well_index] + rest_gas) > E_distribution[well_index]:
                    # finalize the injection to this well
                    inject_to_this_well = E_distribution[well_index] - wells_injected_external_volume[well_index]
                else:
                    # not enough gas to finalize the injection to this well this year
                    inject_to_this_well = rest_gas

                rest_gas -= inject_to_this_well

                # duration (float) of the injection
                duration = min(inject_to_this_well / e_daily, days_in_current_year) # actually inject_to_this_well / e_daily is enough, but sometimes duration = 365.00000006 :(

                # start of the injection
                start = current_t

                # final moment of the injection
                final = duration + start

                # how many days does it take to inject to the current well?
                days_interval = 0

                # if injection was finished in the middle of the day
                if int(final) != final:
                    days_interval += 1

                # account 'full' days
                days_interval += int(final) - int(start)

                # calculate the duration of injection in the first day
                first_day_duration = int(start) + 1 - start

                # calculate the duration of injection in the last day
                if final != int(final):
                    last_day_duration = final - int(final)
                else:
                    last_day_duration = 1.0

                days_duration = [1.0 for d in range(days_interval)]
                days_duration[0] = first_day_duration
                days_duration[-1] = last_day_duration

                # calculate daily injection
                for d in range(days_interval):

                    # find closest point with current E
                    E_closest = min(well.E, key=lambda y: abs(y - wells_injected_external_volume[well_index]))
                    E_ind = (well.E).index(E_closest)
                    Q = well.Q_E[E_ind]
                    I_ind = (well.Q_I).index(Q)

                    if I_ind == 0:
                        dg_rev_dE = (well.I[I_ind+1] - well.I[I_ind]) / (
                            (well.I[I_ind+1] - well.P_I[I_ind+1]) - (well.I[I_ind] - well.P_I[I_ind]))
                        dg_dI = (well.P_I[I_ind+1] - well.P_I[I_ind]) / (well.I[I_ind+1] - well.I[I_ind])
                        dH_dE = (well.Q_I[I_ind+1] - well.Q_I[I_ind]) / (
                            (well.I[I_ind+1] - well.P_I[I_ind+1]) - (well.I[I_ind] - well.P_I[I_ind]))
                    else:
                        dg_rev_dE = (well.I[I_ind] - well.I[I_ind - 1]) / (
                            (well.I[I_ind] - well.P_I[I_ind]) - (well.I[I_ind - 1] - well.P_I[I_ind - 1]))
                        dg_dI = (well.P_I[I_ind] - well.P_I[I_ind-1]) / (well.I[I_ind] - well.I[I_ind-1])
                        dH_dE = (well.Q_I[I_ind] - well.Q_I[I_ind - 1]) / (
                            (well.I[I_ind] - well.P_I[I_ind]) - (well.I[I_ind - 1] - well.P_I[I_ind - 1]))

                    i_daily = dg_rev_dE * e_daily
                    p_daily = dg_dI * i_daily
                    q_daily = dH_dE * e_daily

                    if d == 0 and days_duration[d] != 1.0:
                        # we multiply by days_duration[d] to weight the volume for cases where we change wells during the day
                        q[-1] += q_daily * days_duration[d]
                        p[-1] += p_daily * days_duration[d]
                        i[-1] += i_daily * days_duration[d]
                        e[-1] += e_daily * days_duration[d]
                    else:
                        # we multiply by days_duration[d] to weight the volume for cases where we change wells during the day
                        q.append(q_daily * days_duration[d])
                        p.append(p_daily * days_duration[d])
                        i.append(i_daily * days_duration[d])
                        e.append(e_daily * days_duration[d])

                    wells_injected_external_volume[well_index] += e_daily * days_duration[d]

                if print_log:
                    print('\tInject', inject_to_this_well, ' to well ', wells[well_index].name, ': rest gas = ', rest_gas, ': days_interval = ', days_interval)

                # prepare counter for the next well in this year
                current_t = final

            else: # all gas has been injected to this well
                pass


            # is there still gas left in this year?
            if abs(rest_gas) < volume_error:
                break

        if len(q) != 0:
            for tt in range(len(q)):
                q_field.append(q[tt])
                p_field.append(p[tt])
                i_field.append(i[tt])
                e_field.append(e[tt])

    # make monthly reports and annual volumes
    print('Creating profiles...')

    time = [] # 1,2,3,...7305
    dates = [] # string type: '01.01.2021', '02.01.2021', ...
    mer_dates = [] # string type: '01.2021', '02.2021', ...
    mer_time = [] # 1,2,... 240
    year_dates = [] # 2021, 2022, ...

    q_mer = []
    q_ann = [0.]
    i_mer = []
    i_ann = [0.]
    p_mer = []
    p_ann = [0.]
    e_mer = []
    e_ann = [0.]

    # local summators
    month_duration = 0.
    sum_q = 0.
    sum_i = 0.
    sum_p = 0.
    sum_e = 0.

    # cut MER and annual profiles from cumulative profiles
    current_date = datetime.date(first_year, 1, 1)  # datetime.timedelta(days=1)

    year_dates.append(current_date.year)

    for t in range(len(e_field)):

        time.append(t + 1)
        dates.append(str(current_date.day) + '.' + str(current_date.month) + '.' + str(current_date.year))

        sum_q += q_field[t]
        sum_i += i_field[t]
        sum_p += p_field[t]
        sum_e += e_field[t]
        month_duration += 1

        # if this day is a final day of the month
        if (current_date + datetime.timedelta(days=1)).month != current_date.month:

            # fix current month date
            mer_dates.append(str(current_date.month) + '.' + str(current_date.year))
            mer_time.append(len(mer_time) + 1)

            # fix average rates
            q_mer.append(sum_q / month_duration)
            i_mer.append(sum_i / month_duration)
            p_mer.append(sum_p / month_duration)
            e_mer.append(sum_e / month_duration)

            # add cumulative volumes to the current annual production and injection
            q_ann[-1] += sum_q
            i_ann[-1] += sum_i
            p_ann[-1] += sum_p
            e_ann[-1] += sum_e

            # cut the annual cumulative volumes
            if (current_date + datetime.timedelta(days=1)).year != current_date.year and t != (len(e_field) - 1):
                q_ann.append(0.)
                i_ann.append(0.)
                p_ann.append(0.)
                e_ann.append(0.)
                year_dates.append(current_date.year + 1)

            # reset summators for the next month
            sum_q = 0.
            sum_i = 0.
            sum_p = 0.
            sum_e = 0.
            month_duration = 0.

        # next step
        current_date += datetime.timedelta(days=1)

    # calculate how many compressors do we need using different stepsizes: day, month and year
    one_compr_max = [inj_annual_constraint / volume_unit / 365., inj_annual_constraint / volume_unit / 365., inj_annual_constraint / volume_unit ] # in relative units
    compr_num_float = [max(i_field) / one_compr_max[0], max(i_mer) / one_compr_max[1], max(i_ann) / one_compr_max[2]]
    compr_num = [int(compr_num_float[0]), int(compr_num_float[1]), int(compr_num_float[2])]
    for cn in range(len(compr_num)):
        if compr_num_float[cn] != compr_num[cn]:
            compr_num[cn] += 1

    print('Compressors number = ', compr_num_float, ' ', compr_num)

    if show_graphs:

        plt.plot(time, q_field, color='r')
        plt.plot(time, i_field, color='b')
        plt.plot(time, p_field, color='cyan')
        plt.plot(time, e_field, color='g')
        plt.plot([time[0], time[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[0] + 1):
            plt.plot([time[0], time[-1]], [k * one_compr_max[0], k * one_compr_max[0]], color='maroon', linestyle='-')
        plt.show()

        plt.plot(mer_time, q_mer, color='r', marker='.')
        plt.plot(mer_time, i_mer, color='b', marker='.')
        plt.plot(mer_time, p_mer, color='cyan', marker='.')
        plt.plot(mer_time, e_mer, color='g', marker='.')
        plt.plot([mer_time[0], mer_time[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[1] + 1):
            plt.plot([mer_time[0], mer_time[-1]], [k * one_compr_max[1], k * one_compr_max[1]], color='maroon', linestyle='-')
        plt.show()

        plt.plot(year_dates, q_ann, color='r', marker='.')
        plt.plot(year_dates, i_ann, color='b', marker='.')
        plt.plot(year_dates, p_ann, color='cyan', marker='.')
        plt.plot(year_dates, e_ann, color='g', marker='.')
        plt.plot([year_dates[0], year_dates[-1]], [0, 0], color='k', linestyle='--')
        for k in range(1, compr_num[2] + 1):
            plt.plot([year_dates[0], year_dates[-1]], [k * one_compr_max[2], k * one_compr_max[2]], color='maroon', linestyle='-')
        plt.show()

    # construct the solution
    result = FieldDynamics([time, dates, mer_dates, mer_time, year_dates], [q_field, i_field, p_field, e_field],
                           [q_mer, i_mer, p_mer, e_mer], [q_ann, i_ann, p_ann, e_ann], compr_num_float, compr_num, [0,0])

    print(len(i_field))

    return result

def CheckCumulativeVolumes(profiles, optimal_solution):

    print(profiles.Q_total, optimal_solution.Qmax)
    print(profiles.I_total, sum(optimal_solution.I_wells))
    print(profiles.P_total, sum(optimal_solution.P_wells))
    print(profiles.E_total, optimal_solution.Einj)

def DCPOptimization(n_points, annual_injection, discount, start_point, wells):

    E = [0.]
    for e in annual_injection:
        E.append(E[-1]+e)
    E.remove(E[0])

    # check
    E_max = 0.
    for well in wells:
        E_max += well.E_lin
    if E[-1] > E_max:
        print("Error of E_lin and E_total")
        exit(1)

    if n_points <= 2:
        n_points = 3

    # define array of final points
    final_points = []
    if E[-1] > wells[0].E_lin:
        final_points.append(np.array([wells[0].E_lin, E[-1] - wells[0].E_lin]))
    else:
        final_points.append(np.array([E[-1],0.]))
    if E[-1] > wells[1].E_lin:
        final_points.append(np.array([E[-1] - wells[1].E_lin, wells[1].E_lin]))
    else:
        final_points.append(np.array([0., E[-1]]))

    dr = (final_points[1] - final_points[0]) / (n_points - 1)

    for p in range(1, n_points - 1):
        final_points.insert(p, final_points[0] + dr * p)

    sub_optimal_strategies = []
    for fp in final_points:
        opt_str = BrutForce(n_strategies=3, annual_injection=annual_injection, discount=discount,
              start_point=np.zeros(2), end_point=fp, wells=wells)
        sub_optimal_strategies.append(opt_str)

    DCP = []
    for strategy in sub_optimal_strategies:
        DCP.append(CalculateDCP(strategy, wells, discount))
    print(DCP)

    ShowDCP(DCP, sub_optimal_strategies, annual_injection,
            [final_points[0][0], final_points[0][0]], [0., final_points[0][1]],
            [0., final_points[-1][0]], [final_points[-1][1], final_points[-1][1]], save_option = True)

    return

################################################################################

def BrutForce(n_strategies, annual_injection, discount, start_point, end_point, wells):

    #check that the solution is not trivial (to inject into one pattern)
    for xi in end_point:
        if abs(xi) < 0.01:
            print("The solution is trivial. Exit.\n")
            return

    # make n_strategies = 2k+1, minimal number = 3
    if n_strategies <= 2:
        n_strategies = 3
    if n_strategies%2 == 0:
        n_strategies+=1

    # define middle strategy
    mid = int(n_strategies / 2)
    tg = end_point[1] / end_point[0]

    # list of all strategies (or <=> all trajectories in Eij axes)
    strategies = [ [start_point] for i in range(n_strategies) ]

    # add intermediate points to strategies
    cumulative_total_injection = 0.
    for t in range(len(annual_injection) - 1):
        cumulative_total_injection += annual_injection[t]

        # for middle strategy
        middle_point = np.array([ cumulative_total_injection / (1.+tg), cumulative_total_injection / (1.+tg) * tg ])
        strategies[0].append(middle_point)

        # for lower strategies
        lower_point = np.array([cumulative_total_injection, 0.])
        if cumulative_total_injection > end_point[0]:
            lower_point = np.array([end_point[0], cumulative_total_injection - end_point[0]])
        dr = (lower_point - middle_point) / mid
        for i in range(1, mid+1):
            strategies[i].append(middle_point + i * dr)

        # for upper strategies
        upper_point = np.array([0., cumulative_total_injection])
        if cumulative_total_injection > end_point[1]:
            upper_point = np.array([cumulative_total_injection - end_point[1], end_point[1]])
        dr = (upper_point - middle_point) / mid
        for i in range(1, mid + 1):
            strategies[i+mid].append(middle_point + i * dr)
    for strategy in strategies:
        strategy.append(end_point)
    # now, trajectories for strategies are defined

    # stock_x = []
    # stock_y = []
    # for strategy in strategies:
    #     for point in strategy:
    #         stock_x.append(point[0])
    #         stock_y.append(point[1])
    # plt.scatter(stock_x, stock_y)
    # plt.show()

    # calculate Q_discounted for trajectories
    DCP = []
    for strategy in strategies:
        DCP.append( CalculateDCP(strategy, wells, discount) )

    print(DCP)

    optimum_ind = DCP.index(max(DCP))

    # show DCP
    ShowDCP(DCP, strategies, annual_injection, [end_point[0], end_point[0]], [0, end_point[1]],
            [0, end_point[0]], [end_point[1], end_point[1]], save_option = False)

    return strategies[optimum_ind]

def CalculateDCPandRF(strategy, wells, discount):

    dcp = []
    rf = []

    # calculate total cumulative production for each point of strategy
    for point in strategy:

        local = 0.
        for i,w in enumerate(wells):
            local += w.Q_E_interpolate(float(point[i]))
        rf.append(local)

    # calculate DCP for the trajectory
    for t in range(0, len(strategy)):
        if t == 0:
            dcp.append(0.)
        else:
            d = 1. / (1 + discount) ** (t - 0.5)
            dcp.append(dcp[-1] + d * (rf[t] - rf[t - 1]))

    return dcp, rf

def CalculateDCP(strategy, wells, discount):
    dcp = 0.
    cum_prod = []

    # calculate total (well1 + well2) cumulative production for each point of strategy
    for point in strategy:
        #p = Point([point[i] for i in range(len(wells))], wells)
        local = 0
        for i,w in enumerate(wells):
            local += w.Q_E_interpolate(float(point[i]))
        #cum_prod.append(sum(p.Q))
        cum_prod.append(local)

    # calculate DCP for the trajectory
    for t in range(1, len(strategy)):
        d = 1. / (1 + discount) ** (t - 0.5)
        dcp += d * (cum_prod[t] - cum_prod[t - 1])

    return dcp

def ShowDCP(DCP, strategies, annual_injection, xx1, yy1, xx2, yy2, save_option):
    f, ax = plt.subplots(1)
    min_color = min(DCP)  # 0.297#
    max_color = max(DCP)  # 0.362#
    print(min(DCP), max(DCP))
    cmap = cm.get_cmap('OrRd')  # ('RdYlGn')#('bwr')
    plt.plot(xx1, yy1, 'k--')
    plt.plot(xx2, yy2, 'k--')
    # ax.set_xlim(0, 0.75)
    # ax.set_ylim(0, 0.45)
    ax.set_xticks(np.arange(0, 1., step=0.1))
    ax.set_yticks(np.arange(0, 1., step=0.1))
    E = 0.
    for t in range(len(annual_injection)):
        E += annual_injection[t]
        plt.plot([0, E], [E, 0], 'k-')
    for s, strategy in enumerate(strategies):
        x = []
        y = []
        color = cmap((DCP[s] - min_color) / (max_color - min_color))
        for point in strategy:
            x.append(point[0])
            y.append(point[1])
        plt.plot(x, y, c=color, linewidth=6.0)
    if save_option == True:
        plt.savefig('ShowDCP.pdf', format='pdf')
    plt.show()


################################################################################
def DynProg(E_annual, YearsNum, E_opt, wells, r, discr_E):

    max_edge = math.sqrt(E_opt[0]**2 + E_opt[0]**2)
    dE = max_edge / discr_E
    print("dE = ", dE)



    DP_points = []

    point_last = Point(E_opt, wells)
    DP_points.insert(0, [point_last])
    for i, w in enumerate(wells):
        print("E opt ", str(i + 1), " = ", DP_points[-1][0].E[i], "\tQ opt ", str(i + 1), " = ", DP_points[-1][0].Q[i])


    # go back
    for main_iteration_counter in range(YearsNum - 1, -1, -1):

        #########################
        # 1) create a new points

        # cumulative external injection
        E_current = sum(E_annual[:main_iteration_counter])

        # find point according to E1 + E2 = E_current
        E_edge = [[0.0, 0.0], [0.0, 0.0]]

        for i,w in enumerate(wells):
            if E_current > E_opt[i]:
                E_edge[i][i] = E_opt[i]
                E_edge[i][(i+1)%2] = E_current - E_opt[i]
            else:
                E_edge[i][i] = E_current
                E_edge[i][(i+1)%2] = 0.0

        #print(E_edge)

        delta_edge = math.sqrt((E_edge[0][0] - E_edge[1][0])**2 + (E_edge[0][1] - E_edge[1][1])**2)
        #print("de = ", delta_edge)
        current_points_number = round(delta_edge / dE + 0.5) + 1
        #print("cpn = " , current_points_number)
        if delta_edge != 0.0:
            dE_current = delta_edge / (current_points_number - 1)
        else:
            dE_current = 0.0


        # define Points for current iteration
        current_points = []
        dr = [E_edge[1][0] - E_edge[0][0], E_edge[1][1] - E_edge[0][1]]
        for p in range(current_points_number):
            if current_points_number != 1:
                current_points.append(Point([E_edge[0][0] + p * dr[0] / (current_points_number - 1), E_edge[0][1] + p * dr[1] / (current_points_number - 1)], wells) )
            else:
                current_points.append(Point([E_edge[0][0], E_edge[0][1]], wells))

        DP_points.insert(0, current_points)

        #########################
        # 2) compute optimal value
        for point in DP_points[0]:

            if point.optimal_value != -1:
                print('error in opt value\n')
                exit(1)
            if point.optimal_pointer != -1:
                print('error in opt ptr\n')
                exit(1)

            for j, next_point in enumerate(DP_points[1]):
                if ItCanBeTheNextPoint(point, next_point):
                    value = 1. / (1. + r)**(main_iteration_counter + 0.5) * sum( [ (next_point.Q[i] - point.Q[i]) for i in range(len(wells))] ) + next_point.optimal_value
                    if value > point.optimal_value:
                        point.optimal_value = value
                        point.optimal_pointer = j



    f, ax = plt.subplots(1)
    #ax.set_aspect(1)
    #ax.set_xlim(0, sum(E_annual))
    #ax.set_ylim(0, sum(E_annual))
    plt.xticks(np.arange(0., sum(E_annual), 0.1))
    plt.yticks(np.arange(0., sum(E_annual), 0.1))
    ax.set_xlim(0, 0.75)
    ax.set_ylim(0, 0.4)

    plt.plot([E_opt[0], E_opt[0]], [0, E_opt[1]], 'k--')
    plt.plot([0, E_opt[0]], [E_opt[1], E_opt[1]], 'k--')
    #print(len(DP_points))
    for iteration in DP_points:
        #print("\t", len(iteration))
        for point in iteration:
            #print(point.E)
            plt.plot([point.E[0]], [point.E[1]], color='k', marker='.')

    # plot optimal way
    # go straight
    current_pointer = DP_points[0][0].optimal_pointer
    r0 = [DP_points[0][0].E[0], DP_points[0][0].E[1]]
    # print(current_pointer)
    for main_iteration_counter in range(1, YearsNum+1):
        # print(DP_points[main_iteration_counter][current_pointer].optimal_pointer)
        r1 = [DP_points[main_iteration_counter][current_pointer].E[0],
                DP_points[main_iteration_counter][current_pointer].E[1]]
        plt.plot([r0[0], r1[0]], [r0[1], r1[1]], 'r-')
        #plt.arrow(r0[0], r0[1], r1[0] - r0[0], r1[1] - r0[1], color = 'r', ls='-', lw = 2)
        current_pointer = DP_points[main_iteration_counter][current_pointer].optimal_pointer
        r0 = r1

    #plt.savefig('figure3.pdf', format='pdf')
    #plt.show()
    plt.show()


################################################################################

def ItCanBeTheNextPoint(point, next_point):
    for i in range(len(point.E)):
        if next_point.E[i] < point.E[i]:
            return False

    return True

################################################################################

def ConstructData():

    [I, Q_I, P_I, E, Q_E] = MainDefines(S1, A1, I1, A2, I2)

    MainShow(I, Q_I, P_I, E, Q_E)


################################################################################
def MainDefines(S1, A1, I1, A2, I2):

    # defines
    I = [ind * dI for ind in range(int(I_max / dI))]
    Q_I = f_Q_I(S1, A1, I1, I)
    P_I = f_P_I(A2, I2, I)

    E = [ I[i] - P_I[i] for i in range(len(I))]
    Q_E = copy.copy(Q_I)

    if I_max > E[-1]:
        n_local = (I_max - E[-1]) / dI
        for j in range( int(n_local) + 1):
            E.append(E[-1] +
                     dI)
            Q_E.append(Q_E[-1])
        if E[-1] != I_max:
            E.append(I_max)
            Q_E.append(Q_E[-1])

    return [I, Q_I, P_I, E, Q_E]

################################################################################

def f_Q_I(s1, a1, i1, i):

    result = []
    for ii in i:
        if ii < i1:
            result.append(0.)
        else:
            result.append( s1 * (1. - math.exp(-1. * a1 * (ii - i1)**2)))
    return result

################################################################################

def f_P_I(a2, i2, i):

    result = []
    for ii in i:
        if ii < i2:
            result.append(0.)
        else:
            result.append( 1. / a2 * (math.exp(-1. * a2 * (ii - i2)) - 1) + ii - i2)
    return result

################################################################################

def MainShow(I, Q_I, P_I, E, Q_E):

    label_x_V = "Pore Volume Injected Cumulative, Rpv"
    label_y_Q_I = "Oil Cumulative, Rpv"
    label_y_Q_E = "Oil Cumulative (vs. External Gas Source), Rpv"
    label_y_V = "Gas Cumulative, Rpv"
    label_y_Qr = "Oil Year Rate, Rpv"
#   label_y_Vr = "Gas Year Rate, Rpv"

    host = host_subplot(211, axes_class=AA.Axes)
    par1 = host.twinx()
    host.set_xlabel(label_x_V)
    host.set_ylabel(label_y_Q_I)
    par1.set_ylabel(label_y_V)

    p1, = host.plot(I, Q_I, label=label_y_Q_I)
    p2, = host.plot(E, Q_E, label=label_y_Q_E)
    p3, = par1.plot(I, P_I, label=label_y_V)
    host.legend(loc=4)





    plt.show()

################################################################################

def GetElinear(E, Q_E):

    N = len(E)
    for i in range(N):
        #print(i)
        if abs(Q_E[N - 1 - i] - Q_E[N - 1 - i - 1]) > 0.000001:
            return E[N - 1 - i]

################################################################################

def GetEmax(wells, mult):
    E = 0.
    for w in wells:
        E = E + w.E_lin
    return (E * mult)
################################################################################

################################################################################
################################################################################
if __name__ == '__main__':
    main()
