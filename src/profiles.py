import copy
import datetime

import matplotlib.pyplot as plt

from chatylka_data.chatylka_constants import *

class FieldDynamics:
    def __init__(self, times, data_tr, data_mer, data_ann,
                 compr_num_float, compr_num, max_compressors_number_at_work, exceed_flags, i_mer_wells):

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

        # discounted cumulative production (just for interest)
        self.DCP = sum([self.q_ann[i] / ((1 + r) ** (0.5 + i)) for i in range(len(self.year_dates))])

def CalculateProfiles2(wells_all, sequence, E_distribution_all, e_annual, first_year, print_log,
                       show_graphs, max_compr_injection_annual, max_compressors_number):
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
        raise ValueError('Error in lengths in CalculateProfiles')

    # daily injection limit
    max_i_daily = max_compr_injection_annual * max_compressors_number / volume_unit / 365.

    # how many wells we should inject
    wells_num = 0
    for i in range(len(wells_all)):
        if sequence[i] is not None:
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
