# total volume of external gas in reservoir m3 (RM3) which will have been injected by the end of the project
# it is equal to sum(portion)
volume_unit = 4862682.028

# Gas Formation Volume Factor = reservoir volume / surface volume
# But it's mistake here: I use 1/Bg instead of Bg
Bg = 232 # then real Bg = 1/232

# Oil Formation Volume Factor = reservoir volume / surface volume
Bo = 1.2
# oil density
rho = 0.84

# The number of compressors are used in the oilfiled
Ncompr = 2

# Maximum injectivity produced by one compressor.
# Is it annual, daily, or monthly value???
# Here this value is calculated in reservoir conditions (i.e. / Bg)
max_compr_injection = 135. * 1000. * 1000. / Bg

####################################################################################################
# Data for E(t)
# This is gas supplies in RM3 / Year
portion = [
        827227.5639,
        733917.6983,
        738013.3463,
        524818.9293,
        388085.5369,
        300185.7802,
        237263.0647,
        192122.7119,
        161535.0881,
        133575.1259,
        112187.3892,
        95257.64978,
        85030.34537,
        72453.9945,
        62022.90808,
        53437.12403,
        45683.95528,
        38945.20657,
        33138.48287,
        27780.12683
    ]

# isn't used here (special case where gas supplies are equal every year)
portion_aver = [
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014,
        243134.1014
]

# Special multiplier for gas supplies scaling
E_mult = 1.0
#YearsNum = len(portion)


####################################################################################################
# Data for WELLS

# Directory where all txt-files Q_I, P_I and Q_E are there
data_dir = './'

# Postfixes for file-names (each postfix means the real well's name)
names = [
        '10',
        '11',
        '12',
        '13',
        '15',
        '17',
        '18',
        '20',
        '24',
        '31',
        '40',
        '45',
        '46',
        '51'
    ]

# Mask for well names (mask can be used for articles where real well names can't be used)
names_mask = {
        '10': '#1',
        '11': '#2',
        '12': '#3',
        '13': '#4',
        '15': '#5',
        '17': '#6',
        '18': '#7',
        '20': '#8',
        '24': '#9',
        '31': '#10',
        '40': '#11',
        '45': '#12',
        '46': '#13',
        '51': '#14'
    }

# How wells are named in reservoir simulation model
gdm_mask = {
        '10': '10ST2',
        '11': '11ST3',
        '12': '12',
        '13': '13',
        '15': '15',
        '17': '17',
        '18': '18',
        '20': '20',
        '24': '24',
        '31': '31',
        '40': '40ST2',
        '45': '45',
        '46': '46',
        '51': '51'
}

# Current injectivity for the injection wells
WCONINJH = {
        '10': 72.767,
        '11': 185.5,
        '12': 162.6,
        '13': 277.615,
        '15': 330.259,
        '17': 246.7,
        '18': 265.607,
        '20': 241.033,
        '24': 329.333,
        '31': 447.071,
        '40': 194.8,
        '45': 109.933,
        '46': 376.4,
        '51': 189.567
        }


# Current well efficiencies for the injection wells
WEFAC = {
        '10': 1.,
        '11': 1.,
        '12': 1.,
        '13': 0.161,
        '15': 0.226,
        '17': 0.065,
        '18': 0.323,
        '20': 0.387,
        '24': 0.065,
        '31': 0.433,
        '40': 0.6,
        '45': 0.554,
        '46': 0.133,
        '51': 0.901
}

##################################
# Economics

netback = 13867 # rub/t
netback_ndd = 19895 # rub/t
royalty = 4631 # rub/t # or 7748
royalty_ndd = 4631 # 2740 # rub/t
income_tax = 0.20 # 20%
ndd_coeff = 0.50 # 50%
c_oil = 74 # rub/t
c_gas = 0.01 # rub/m3
niz_lu_ndd = 1075000 # t


first_year = 2021

# discount rate (annual)
r = 0.14

# The year the economics is reduced to.
first_year_econ = 2018 # period = 0.5

