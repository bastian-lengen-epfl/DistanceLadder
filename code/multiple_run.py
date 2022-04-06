'''
This script takes the pre-processed data from the ./data_tmp/ folder, fit them according to the fit_parameters.py file,
print,  save the result of the fit and draw a few plots.
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
from relativistic_corrections import RLB_correction, interpolated_K_corr_Cep, interpolated_K_corr_TRGB
from outliers_rejection import single_kappa_clipping
from plots import plot_individual_PL, plot_global_PL, plot_SNe


def single_K_corr_run(fit_dir, **kwargs):
    ### Create the fit directory
    fit_dir = fit_dir + '/'
    if not os.path.exists(fit_dir):
        print(f'I will create the {fit_dir} directory for you !')
        os.mkdir(fit_dir)

    # Load the data
    DF_dict, SNe_other, table_Kcorr_Cep, table_Kcorr_TRGB = load_data(data_tmp='./data_tmp/',
                                                                      data_static='./data_static/')
    interpolated_fun_list = kwargs.get('interpolated_fun_list')
    interpolated_fun = kwargs.get('interpolated_fun')

    ### Relativistic corrections
    # Cepheids relativistic corrections
    if fp.include_Cepheids == True:
        if fp.RLB_correction == True:
            print('Correcting Cepheids for the RLB...')
            RLB_correction(DF_dict)
        if fp.Kcorr_Cep == True:
            print('Start of the K-correction for the Cepheids...')
            interpolated_fun_list = interpolated_K_corr_Cep(DF_dict, table_Kcorr_Cep, kwargs.get('EBV_Cep', fp.EBV_Cep),
                                                            interpolated_fun_list)
    # TRGB relativistic corrections
    if fp.include_TRGB == True:
        if fp.Kcorr_TRGB == True:
            print('Start of the K-correction the TRGB...')
            interpolated_fun = interpolated_K_corr_TRGB(DF_dict, table_Kcorr_TRGB,
                                                        kwargs.get('Teff_TRGB', fp.Teff_TRGB),
                                                        kwargs.get('logg_TRGB', fp.logg_TRGB),
                                                        kwargs.get('FeH_TRGB', fp.FeH_TRGB),
                                                        kwargs.get('EBV_TRGB', fp.EBV_TRGB),
                                                        interpolated_fun)

    ### Fitting
    # Outliers rejection (Cepheids and SNe only)
    if (fp.include_Cepheids == True) or (fp.fit_aB == True):
        if fp.outlier_rejection == True:
            print('Running the outlier rejection...')
            if fp.fit_aB == False:
                SNe_other = pd.DataFrame()  # Emtpy SNe_other if no SNe
            DF_dict, DF_dict_outliers = single_kappa_clipping(DF_dict, SNe_other, kappa=2.7, work_dir=fit_dir)
    # Fit
    print('Computing H0...')
    y, q_dict, L, Sigma2 = fit_distance_ladder(DF_dict)

    ### Save setup
    print('Saving the fit_parameters.py used...')
    shutil.copyfile('./code/fit_parameters.py', fit_dir + 'fit_parameters.py')
    print('Saving the datasets used...')
    data_dir = fit_dir + 'data/'
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
    y_pkl = open(fit_dir + 'y.pickle', 'wb')
    pickle.dump(y, y_pkl)
    y_pkl.close()
    q_dict_pkl = open(fit_dir + 'q_dict.pickle', 'wb')
    pickle.dump(q_dict, q_dict_pkl)
    q_dict_pkl.close()
    L_pkl = open(fit_dir + 'L.pickle', 'wb')
    pickle.dump(L, L_pkl)
    L_pkl.close()
    Sigma2_pkl = open(fit_dir + 'Sigma2.pickle', 'wb')
    pickle.dump(Sigma2, Sigma2_pkl)
    Sigma2_pkl.close()

    ### Plots
    if fp.show_plots == True:
        if fp.outlier_rejection == False:
            if fp.include_Cepheids == True:
                plot_individual_PL(DF_dict, q_dict, dict({}), work_dir=fit_dir)  # Empty dict for outliers
                plot_global_PL(DF_dict, q_dict, dict({}), work_dir=fit_dir)  # Empty dict for outliers
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, dict({}), work_dir=fit_dir)
        else:
            if fp.include_Cepheids == True:
                plot_individual_PL(DF_dict, q_dict, DF_dict_outliers, work_dir=fit_dir)
                plot_global_PL(DF_dict, q_dict, DF_dict_outliers, work_dir=fit_dir)
            if fp.fit_aB == True:
                plot_SNe(DF_dict, q_dict, DF_dict_outliers, work_dir=fit_dir)

    ### Display results
    for str in q_dict:
        print(f'{str} : {q_dict[str]}')

    return q_dict, interpolated_fun_list, interpolated_fun

def multiple_run(fit_name, work_dir):
    ### Check for inconsiticency in fit_parameters
    if (fp.include_Cepheids == False) and (fp.multiple_Cep == True):
        print('ERROR: Cannot run multiple Cepheids with include_Cepheids=False.')
        return 1
    if (fp.include_TRGB == False) and (fp.multiple_TRGB == True):
        print('ERROR: Cannot run multiple TRGB with include_TRGB=False.')
        return 1
    if (fp.multiple_Cep == True) and (fp.multiple_TRGB == True):
        print('ERROR: Cannot run both multiple Cepheids and multiple TRGB at the same time.')
        return 1
    if (fp.multiple_Cep == False) and (fp.multiple_TRGB == False):
        print('ERROR: At least one of multiple_Cep or multiple_TRGB has to be True.')
        return 1
    if (fp.multiple_Cep == True) and (fp.Kcorr_Cep == False):
        print('ERROR: Cannot run multiple Cepheids with Kcorr_Cep = False.')
        return 1
    if (fp.multiple_TRGB == True) and (fp.Kcorr_TRGB == False):
        print('ERROR: Cannot run multiple TRGB with Kcorr_TRGB = False.')
        return 1

    ### Create the main working directory
    if not os.path.exists(work_dir):
        print(f'I will create the {work_dir} directory for you !')
        os.mkdir(work_dir)
    work_dir = work_dir + fit_name + '/'
    if not os.path.exists(work_dir):
        print(f'I will create the {work_dir} directory for you !')
        os.mkdir(work_dir)



    ### Multiple K-corr run
    interpolated_fun_list, interpolated_fun = None, None
    if fp.multiple_Cep == True:
        q_summary = []
        for EBV in fp.EBV_Cep_multi:
            print(f'*'.center(80, '*'))
            print(f' EBV = {EBV} '.center(80, '*'))
            print(f'*'.center(80, '*'))
            q_string, interpolated_fun_list, interpolated_fun = \
                single_K_corr_run(work_dir+f'EBV_Cep={EBV}', EBV_Cep=EBV, interpolated_fun_list=interpolated_fun_list)
            q_string['EBV_Cep'] = EBV
            q_summary.append(q_string)
    elif fp.multiple_TRGB == True:
        q_summary = []
        for Teff in fp.Teff_TRGB_multi:
            for EBV in fp.EBV_TRGB_multi:
                print(f'*'.center(80, '*'))
                print(f' Teff = {Teff}, EBV = {EBV} '.center(80, '*'))
                print(f'*'.center(80, '*'))
                q_string, interpolated_fun_list, interpolated_fun = \
                    single_K_corr_run(work_dir+f'Teff_TRGB={Teff}_EBV_TRGB={EBV}', Teff_TRGB=Teff, EBV_TRGB=EBV,
                                      logg_TRGB=fp.logg_TRGB, FeH_TRGB=fp.FeH_TRGB, interpolated_fun=interpolated_fun)
                q_string['Teff_TRGB'] = Teff
                q_string['EBV_TRGB']  = EBV
                q_string['logg_TRGB'] = fp.logg_TRGB
                q_string['FeH_TRGB'] = fp.FeH_TRGB
                q_summary.append(q_string)


    ### Save pkl
    print('Pickling the summary variable q_summary...')
    q_summary_pkl = open(work_dir + 'q_summary.pickle', 'wb')
    pickle.dump(q_summary, q_summary_pkl)
    q_summary_pkl.close()
    return 0

if __name__=='__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description='TO DO',
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument(type=str, dest='fit_name', metavar='fit_name',
                        help='Give a name to your fit. Determines the name of the folder containing all the results'
                             'from the fit.')
    parser.add_argument('--dir', type=str, dest='work_dir', metavar='work_dir',
                        default='./work_dir/',
                        help='Name of the working directory')

    args = parser.parse_args()
    multiple_run(args.fit_name, args.work_dir)

