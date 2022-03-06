"""
This module contains a function that loads the previously pre-processed data and return the corresponding DataFrames
that are used by the other functions.
"""
import pandas as pd
import fit_parameters as fp

def load_data(work_dir='./'):
    '''
    Loads and return DF_dict, a dictionary of pandas DataFrame of the different .csv. It also return the SNe that are
    out of the [z_min, z_max] interval in a SNe_other DataFrame

    :type   work_dir: string
    :param  work_dir: working directory, by default ./
    '''
    Data_dir = work_dir + 'data/'

    # Create an empty dictionnary of dataset
    DF_dict={}

    # Loads the Cepheids related DataFrame
    if fp.include_Cepheids==True:
        print('Loading Cepheids & Cepheids anchors...')
        DF_dict['Cepheids'] = pd.read_csv(Data_dir+'Cepheids.csv', sep=',')
        DF_dict['Cepheids'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I']
        DF_dict['Cepheids_anchors'] = pd.read_csv(Data_dir + 'Cepheids_anchors.csv', sep=',')
        DF_dict['Cepheids_anchors'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'mu', 'sig_mu']
        if fp.include_MW==True:
            print('Loading MW Cepheids...')
            DF_dict['Cepheids_MW'] = pd.read_csv(Data_dir+'Cepheids_MW.csv', sep=',')
            DF_dict['Cepheids_MW'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'pi', 'sig_pi']
        print('Loading SNe for the Cepheids-host galaxies...')
        DF_dict['SNe_Cepheids'] = pd.read_csv(Data_dir + 'SNe_Cepheids.csv', sep=',')
        DF_dict['SNe_Cepheids'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the TRGB related DataFrame
    if fp.include_TRGB==True:
        print('Loading TRGB & TRGB anchors...')
        DF_dict['TRGB'] = pd.read_csv(Data_dir+'TRGB.csv', sep=',')
        DF_dict['TRGB'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z']
        DF_dict['TRGB_anchors'] = pd.read_csv(Data_dir + 'TRGB_anchors.csv', sep=',')
        DF_dict['TRGB_anchors'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z', 'mu', 'sig_mu']
        print('Loading SNe for the TRGB-host galaxies...')
        DF_dict['SNe_TRGB'] = pd.read_csv(Data_dir + 'SNe_TRGB.csv', sep=',')
        DF_dict['SNe_TRGB'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the SNe for the aB fit
    if fp.fit_aB==True:
        print('Loading SNe for the aB fit...')
        DF_dict['SNe_Hubble'] = pd.read_csv(Data_dir + 'SNe_Hubble.csv', sep=',')
        DF_dict['SNe_Hubble'].columns = ['name', 'mB', 'sig_mB', 'z', 'sig_z']
        filter = ((DF_dict['SNe_Hubble']['z'] > fp.z_min) & (DF_dict['SNe_Hubble']['z'] < fp.z_max))
        SNe_other = DF_dict['SNe_Hubble'][~filter]
        DF_dict['SNe_Hubble'] = DF_dict['SNe_Hubble'][filter].reset_index(drop=True)

    return DF_dict, SNe_other
