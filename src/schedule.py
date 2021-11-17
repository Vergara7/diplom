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
