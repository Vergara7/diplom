class Report:
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
