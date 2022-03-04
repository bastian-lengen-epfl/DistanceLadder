"""
This module contains a functions that loads the previously pre-processed data and return the corresponding dataset
that is used by the other functions.
"""
import pandas as pd
import Fit_parameters as Fp

def load_data(work_dir='./'):
    '''
    Loads and return a dictionary of pandas DataFrame of the different .csv

    :type   work_dir: string
    :param  work_dir:
    '''
    Data_dir = work_dir + 'Data/'

    # Create an empty dictionnary of dataset
    DF_dict={}

    # Loads the Cepheids related DataFrame
    if Fp.include_Cepheids==True:
        print('Loading Cepheids & Cepheids anchors...')
        DF_dict['Cepheids'] = pd.read_csv(Data_dir+'Cepheids.csv', sep=',')
        DF_dict['Cepheids'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I']
        DF_dict['Cepheids_anchors'] = pd.read_csv(Data_dir + 'Cepheids_anchors.csv', sep=',')
        DF_dict['Cepheids_anchors'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'mu', 'sig_mu']
        if Fp.include_MW==True:
            print('Loading MW Cepheids...')
            DF_dict['Cepheids_MW'] = pd.read_csv(Data_dir+'Cepheids_MW.csv', sep=',')
            DF_dict['Cepheids_MW'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'pi', 'sig_pi']
        print('Loading SNe for the Cepheids-host galaxies...')
        DF_dict['SNe_Cepheids'] = pd.read_csv(Data_dir + 'SNe_Cepheids.csv', sep=',')
        DF_dict['SNe_Cepheids'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the TRGB related DataFrame
    if Fp.include_TRGB==True:
        print('Loading TRGB & TRGB anchors...')
        DF_dict['TRGB'] = pd.read_csv(Data_dir+'TRGB.csv', sep=',')
        DF_dict['TRGB'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z']
        DF_dict['TRGB_anchors'] = pd.read_csv(Data_dir + 'TRGB_anchors.csv', sep=',')
        DF_dict['TRGB_anchors'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z', 'mu', 'sig_mu']
        print('Loading SNe for the TRGB-host galaxies...')
        DF_dict['SNe_TRGB'] = pd.read_csv(Data_dir + 'SNe_TRGB.csv', sep=',')
        DF_dict['SNe_TRGB'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the SNe for the aB fit
    if Fp.fit_aB==True:
        print('Loading SNe for the aB fit...')
        DF_dict['SNe_Hubble'] = pd.read_csv(Data_dir + 'SNe_Hubble.csv', sep=',')
        DF_dict['SNe_Hubble'].columns = ['name', 'mB', 'sig_mB', 'z', 'sig_z']

    return DF_dict
