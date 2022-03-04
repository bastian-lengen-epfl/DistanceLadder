import sys
import Fit_parameters as Fp
from Fitting import fit_distance_ladder
from Load_data import load_data
from Relativistic_corrections import RLB_correction, K_corr_Cep, K_corr_TRGB


def main() -> int:
    ### Load the data
    DF_dict = load_data()

    ### Relativistic corrections
    # Cepheids relativistic corrections
    if Fp.include_Cepheids == True:
        if Fp.RLB_correction == True:
            print('Correcting Cepheids for the RLB...')
            DF_dict = RLB_correction(DF_dict)
        if Fp.Kcorr_Cep == True:
            print('K-correcting the Cepheids...')
            DF_dict = K_corr_Cep(DF_dict)
    # TRGB relativistic corrections
    if Fp.include_TRGB == True:
        if Fp.Kcorr_TRGB == True:
            print('K-correcting the TRGB...')
            DF_dict = K_corr_TRGB(DF_dict)

    ### Fitting
    # Outliers rejection (Cepheids and SNe only)
    if (Fp.include_Cepheids == True) or (Fp.fit_aB == True):
        if Fp.outlier_rejection == True:
            pass
    # Fit
    y, q_dict, L= fit_distance_ladder(DF_dict)


    ### Display and save results
    # q results
    for str in q_dict:
        print(f'{str} : {q_dict[str]}')


    return 0

if __name__ == '__main__':
    sys.exit(main())  # next section explains the use of sys.exit