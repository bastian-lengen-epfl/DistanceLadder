'''
This script takes the pre-processed data from the ./data/ folder, fit them according to the fit_parameters.py file,
print and save the results of the fit and finally draw a few plots.
'''
import sys
import os
import argparse as ap
import fit_parameters as fp
import pandas as pd
import shutil
import pickle
from fitting import fit_distance_ladder
from load_data import load_data
from relativistic_corrections import RLB_correction, interpolated_K_corr_Cep,  interpolated_K_corr_TRGB
from outliers_rejection import single_kappa_clipping
from plots import plot_individual_PL, plot_global_PL, plot_SNe


def run(fit_name, work_dir):
    ### Create the working directory:
    if not os.path.exists(work_dir):
        print(f'I will create the {work_dir} directory for you !')
        os.mkdir(work_dir)
    work_dir = work_dir + fit_name + '/'
    if not os.path.exists(work_dir):
        print(f'I will create the {work_dir} directory for you !')
        os.mkdir(work_dir)

    ### Load the data
    DF_dict, SNe_other, table_Kcorr = load_data(data='./data/', data_static='./data_static/')

    ### Relativistic corrections
    # Cepheids relativistic corrections
    if fp.include_Cepheids == True:
        if fp.RLB_correction == True:
            print('Correcting Cepheids for the RLB...')
            RLB_correction(DF_dict)
        if fp.Kcorr_Cep == True:
            print('Start of the K-correction for the Cepheids...')
            interpolated_K_corr_Cep(DF_dict, table_Kcorr, EBV=fp.EBV_Cep)
    # TRGB relativistic corrections
    if fp.include_TRGB == True:
        if fp.Kcorr_TRGB == True:
            print('Start of the K-correction the TRGB...')
            interpolated_K_corr_TRGB(DF_dict, table_Kcorr, fp.Teff_TRGB,
                                     fp.logg_TRGB, fp.FeH_TRGB, fp.EBV_TRGB)

    ### Fitting
    # Outliers rejection (Cepheids and SNe only)
    if (fp.include_Cepheids == True) or (fp.fit_aB == True):
        if fp.outlier_rejection == True:
            print('Running the outlier rejection...')
            if fp.fit_aB == False:
                SNe_other = pd.DataFrame() # Emtpy SNe_other if no SNe
            DF_dict, DF_dict_outliers = single_kappa_clipping(DF_dict, SNe_other, kappa=fp.kappa, work_dir=work_dir)
    # Fit
    print('Computing H0...')
    y, q_dict, L, Sigma2 = fit_distance_ladder(DF_dict)

    ### Save setup
    print('Saving the fit_parameters.py used...')
    shutil.copyfile('./code/fit_parameters.py',work_dir+'fit_parameters.py')
    print('Saving the datasets used...')
    data_dir = work_dir + 'data/'
    if not os.path.exists(data_dir):
        print(f'I will create the {data_dir} directory for you !')
        os.mkdir(data_dir)
    if fp.include_Cepheids == True:
        DF_dict['Cepheids'].to_csv(data_dir + 'Cepheids.csv', index=False)
        DF_dict['Cepheids_anchors'].to_csv(data_dir + 'Cepheids_anchors.csv', index=False)
        if fp.include_MW == True:
            DF_dict['Cepheids_MW'].to_csv(data_dir + 'Cepheids_MW.csv', index=False)
        DF_dict['SNe_Cepheids'].to_csv(data_dir + 'SNe_Cepheids.csv', index=False)
    if fp.include_TRGB == True:
        DF_dict['TRGB'].to_csv(data_dir + 'TRGB.csv', index=False)
        DF_dict['TRGB_anchors'].to_csv(data_dir + 'TRGB_anchors.csv', index=False)
        DF_dict['SNe_TRGB'].to_csv(data_dir + 'SNe_TRGB.csv', index=False)
    if fp.fit_aB == True:
        DF_dict['SNe_Hubble'].to_csv(data_dir + 'SNe_Hubble.csv', index=False)

    ### Save results
    print('Pickling the variable y, q_dict, L and Sigma2...')
    y_pkl = open(work_dir + 'y.pickle', 'wb')
    pickle.dump(y, y_pkl)
    y_pkl.close()
    q_dict_pkl = open(work_dir + 'q_dict.pickle', 'wb')
    pickle.dump(q_dict, q_dict_pkl)
    q_dict_pkl.close()
    L_pkl = open(work_dir + 'L.pickle', 'wb')
    pickle.dump(L, L_pkl)
    L_pkl.close()
    Sigma2_pkl = open(work_dir + 'Sigma2.pickle', 'wb')
    pickle.dump(Sigma2,Sigma2_pkl)
    Sigma2_pkl.close()

    ### Plots
    if fp.show_plots == True:
        if fp.outlier_rejection == False:
            if fp.include_Cepheids == True:
                plot_individual_PL(DF_dict, q_dict, dict({}), work_dir=work_dir) # Empty dict for outliers
                plot_global_PL(DF_dict, q_dict, dict({}), work_dir=work_dir)  # Empty dict for outliers
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, dict({}), work_dir=work_dir)
        else:
            if fp.include_Cepheids == True:
                plot_individual_PL(DF_dict, q_dict, DF_dict_outliers, work_dir=work_dir)
                plot_global_PL(DF_dict, q_dict, DF_dict_outliers, work_dir=work_dir)
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, DF_dict_outliers, work_dir=work_dir)

    ### Display results
    for str in q_dict:
        print(f'{str} : {q_dict[str]}')

    return 0

if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description='Fits the distance ladder according to the data from the ./data folder '
                                           'and the fit parameters from the ./code/fit_parameters.py file. The results '
                                           'are then saved in the working directory.',
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument(type=str, dest='fit_name', metavar='fit_name',
                        help='Give a name to your fit. Determines the name of the folder containing all the results'
                             'from the fit.')
    parser.add_argument('--dir', type=str, dest='work_dir', metavar='work_dir',
                        default='./work_dir/',
                        help='Name of the working directory.')

    args = parser.parse_args()
    run(args.fit_name, args.work_dir)

