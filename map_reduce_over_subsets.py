import argparse
import itertools
import os

from scipy.optimize import minimize

import yt.wrapper as yt

from pycore.yt_utils import get_yt_client
from pycore.porto_utils import PORTO_LAYER_UBUNTU

data_dir = os.path.join(os.getcwd(), 'diplom/chatylka_data')
from chatylka_data.chatylka_constants import *

from src.common import *

from src.well import GetWellsFromPreparedTxtFiles
from src.profiles import CalculateProfiles2

def parse_args():
    parser = argparse.ArgumentParser(description='Sequential strategies runner')

    parser.add_argument('-i', '--input-table', required=True, help='input table')
    parser.add_argument('-o', '--output-table', required=True, help='output table')

    return parser.parse_args()

class SequentialStrategiesRunner(object):
    def __init__(self, wells, e_annual, optimize_e_distr=False):
        self._wells = wells
        self._e_annual = e_annual
        self._optimize_e_distr = optimize_e_distr


    def _get_e_distr(self, permutaion):
        used_wells = sum([1 for v in permutaion if v is not None])
        init_distr = [0 if v is None else 1 / used_wells for v in permutaion]
        if not self._optimize_e_distr:
            return init_distr

        def calculateDCP(e_distr):
            distr = [0 for w in self._wells]
            cur = 0
            for i, v in enumerate(permutaion):
                if v is not None:
                    distr[i] = e_distr[cur]
                    cur += 1
            profiles = CalculateProfiles2(
                self._wells, full_perm, distr, self._e_annual, first_year,
                print_log=False, show_graphs=False,
                max_compr_injection_annual=max_compr_injection,
                max_compressors_number=Ncompr)
            profiles.CalculateEconomics(netback, royalty, income_tax, c_oil, c_gas, r, first_year_econ)
            return profiles.DCP

        eq_cons = {
            'type': 'eq',
            'fun': lambda x: sum(x) - 1,
            'jac': lambda x: np.array([1.0 for w in range(self._wells)])
        }

        bounds =

        res = minimize(
            calculateDCP, x0, method='SLSQP', constraints=eq_cons,
            options={'ftol': 1e-8}, bounds=bounds, args=wells)






    def __call__(self, row):
        for permutation in itertools.permutations(row['ids']):
            full_perm = [None for w in self._wells]
            for i, v in enumerate(permutation):
                full_perm[v] = i
            e_distr = self._get_e_distr(full_perm)
            profiles = CalculateProfiles2(
                self._wells, full_perm, e_distr, self._e_annual, first_year,
                print_log=False, show_graphs=False,
                max_compr_injection_annual=max_compr_injection,
                max_compressors_number=Ncompr)
            profiles.CalculateEconomics(netback, royalty, income_tax, c_oil, c_gas, r, first_year_econ)
            yield {'permutaion': permutation, 'NPV': profiles.NPV / 10**9, 'DCP': profiles.DCP}


def main():
    args = parse_args()
    E_annual = [E_mult * p / volume_unit for p in portion]
    E_total = sum(E_annual)
    wells = GetWellsFromPreparedTxtFiles(
        data_dir=data_dir, Q_E_str='Q_E_', Q_I_str='Q_I_', P_I_str='P_I_',
        well_names=names, well_masks=names_mask, init_option ='files',
        E_total=E_total, volume_unit=volume_unit)
    wells.sort(key=lambda well: max(well.Q_E), reverse=True)
    for well in wells:
        well.FindPointOfEffectiviness()

    yt_client = get_yt_client()
    runner = SequentialStrategiesRunner(wells, E_annual)
    print(list(runner({'ids': [0, 1, 2]})))
    return
    yt_client.run_map(runner, args.input_table, args.output_table)
    yt_client.sort(args.input_table, args.output_table)


if __name__ == '__main__':
    main()
