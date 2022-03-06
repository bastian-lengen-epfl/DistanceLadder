'''
This script take the pre-processed data from the ./data folder, fit them according to the fit_parameters.py file,
and print the result for the different parameters.
'''
import sys
import fit_parameters as fp
import pandas as pd
from fitting import fit_distance_ladder
from load_data import load_data
from relativistic_corrections import RLB_correction, K_corr_Cep, K_corr_TRGB
from outliers_rejection import single_kappa_clipping
from plots import plot_individual_PL, plot_global_PL, plot_SNe


def main() -> int:
    ### Load the data
    DF_dict, SNe_other = load_data()

    ### Relativistic corrections
    # Cepheids relativistic corrections
    if fp.include_Cepheids == True:
        if fp.RLB_correction == True:
            print('Correcting Cepheids for the RLB...')
            DF_dict = RLB_correction(DF_dict)
        if fp.Kcorr_Cep == True:
            print('K-correcting the Cepheids...')
            DF_dict = K_corr_Cep(DF_dict)
    # TRGB relativistic corrections
    if fp.include_TRGB == True:
        if fp.Kcorr_TRGB == True:
            print('K-correcting the TRGB...')
            DF_dict = K_corr_TRGB(DF_dict)

    ### Fitting
    # Outliers rejection (Cepheids and SNe only)
    if (fp.include_Cepheids == True) or (fp.fit_aB == True):
        if fp.outlier_rejection == True:
            print('Running the outlier rejection...')
            if fp.fit_aB == False:
                SNe_other = pd.DataFrame() # Emtpy SNe_other if no SNe
            DF_dict, DF_dict_outliers = single_kappa_clipping(DF_dict, SNe_other, kappa=2.7, work_dir='./')
    # Fit
    print('Computing H0...')
    y, q_dict, L= fit_distance_ladder(DF_dict)


    ### Display results
    # q results
    for str in q_dict:
        print(f'{str} : {q_dict[str]}')

    ### Plots
    if fp.show_plots == True:
        if fp.outlier_rejection == False:
            if fp.include_Cepheids == True:
                pass
                plot_individual_PL(DF_dict, q_dict, dict({})) # Empty dict for outliers
                plot_global_PL(DF_dict, q_dict, dict({}))  # Empty dict for outliers
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, dict({}))
        else:
            if fp.include_Cepheids == True:
                pass
                plot_individual_PL(DF_dict, q_dict, DF_dict_outliers)
                plot_global_PL(DF_dict, q_dict, DF_dict_outliers)
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, DF_dict_outliers)
    return 0

if __name__ == '__main__':
    sys.exit(main())

