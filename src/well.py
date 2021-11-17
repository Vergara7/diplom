import os

from scipy.interpolate import interp1d

# I
dI = 0.01
I_max = 10

class Well:

    def __init__(self, name, option, data):
        self.name = name

        assert option == 'files'

        self.I = data[0]
        self.Q_I = data[1]
        self.P_I = data[2]
        self.E = data[3]
        self.Q_E = data[4]

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


def GetElinear(E, Q_E):
    N = len(E)
    for i in range(N):
        if abs(Q_E[N - 1 - i] - Q_E[N - 1 - i - 1]) > 0.000001:
            return E[N - 1 - i]


def MainDefines(S1, A1, I1, A2, I2):
    I = [ind * dI for ind in range(int(I_max / dI))]
    Q_I = f_Q_I(S1, A1, I1, I)
    P_I = f_P_I(A2, I2, I)

    E = [ I[i] - P_I[i] for i in range(len(I))]
    Q_E = copy.copy(Q_I)

    if I_max > E[-1]:
        n_local = (I_max - E[-1]) / dI
        for j in range( int(n_local) + 1):
            E.append(E[-1] + dI)
            Q_E.append(Q_E[-1])
        if E[-1] != I_max:
            E.append(I_max)
            Q_E.append(Q_E[-1])

    return [I, Q_I, P_I, E, Q_E]


def f_Q_I(s1, a1, i1, i):

    result = []
    for ii in i:
        if ii < i1:
            result.append(0.)
        else:
            result.append( s1 * (1. - math.exp(-1. * a1 * (ii - i1)**2)))
    return result


def f_P_I(a2, i2, i):

    result = []
    for ii in i:
        if ii < i2:
            result.append(0.)
        else:
            result.append( 1. / a2 * (math.exp(-1. * a2 * (ii - i2)) - 1) + ii - i2)
    return result

def GetWellsFromPreparedTxtFiles(data_dir, Q_E_str, Q_I_str, P_I_str,
                                 well_names, well_masks, init_option, E_total, volume_unit):
    wells = []
    for i in range(len(well_names)):
        with open(os.path.join(data_dir, Q_E_str + str(well_names[i]) + '.txt'), 'r') as f:
            lines = f.readlines()
        E = []
        Q = []
        for line in lines:
            EQ = line.split('\t')
            E.append(eval(EQ[0]) / volume_unit)
            Q.append(eval(EQ[1]) / volume_unit)

        with open(os.path.join(data_dir, Q_I_str + str(well_names[i]) + '.txt'), 'r') as f:
            lines = f.readlines()
        I = []
        Q_I = []
        for line in lines:
            IQ = line.split('\t')
            I.append(eval(IQ[0]) / volume_unit)
            Q_I.append(eval(IQ[1]) / volume_unit)

        with open(os.path.join(data_dir, Q_I_str + str(well_names[i]) + '.txt'), 'r') as f:
            lines = f.readlines()
        I_ = []
        P_I = []
        for line in lines:
            IP = line.split('\t')
            I_.append(eval(IP[0]) / volume_unit)
            P_I.append(eval(IP[1]) / volume_unit)

        if I != I_:
            print(f'Error in file loading for well {str(well_names[i])}\n')
            exit(1)

        wells.append(Well(well_masks[str(well_names[i])], init_option, [I, Q_I, P_I, E, Q]))

    # check the loaded data
    E_lin_summ = 0.
    for i, w in enumerate(wells):
        E_lin_summ += w.E_lin
    if E_total >= E_lin_summ:
        print('Error: trivial solution: E = E_lin')
        print(f'E_total = {E_total}; but the summ of E_lin = {E_lin_summ}')
        print('Exit...')
        exit(1)
    return wells
