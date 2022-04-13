'''
This module contains a function that loads the previously pre-processed data and return the corresponding DataFrames
that are used by the other functions.
'''
import pandas as pd
import fit_parameters as fp

def load_data(data_tmp='./data_tmp/', data_static='./data_static/') -> bool:
    '''
    Load the data from the data_tmp and data_static folders and put them in different DataFrames that will be later used by the other
    functions.

    Parameters
    ----------
    data_tmp : str
        data_tmp directory, by default ./data_tmp/
    data_static : str
        data_static directory, by default ./data_static/

    Returns
    -------
    DF_dict : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have to be fitted.
    SNe_other : pandas DataFrame
        DataFrame that contains the SNe that are outside the fitting range [z_min, z_max] for the Hubble diagram.
    table_Kcorr : pandas DataFrame
        DataFrame that contains Table 1a from Anderson (2022) [2022A&A...658A.148A] that will be used for the
        K-corrections. Note that it can return None if the fit doesn't apply K-correction (see fit_parameters.py).
    '''
    # Create an empty dictionnary of dataset
    DF_dict=dict()

    # Loads the Cepheids related DataFrame
    if fp.include_Cepheids==True:
        print('Loading Cepheids & Cepheids anchors...')
        DF_dict['Cepheids'] = pd.read_csv(data_tmp+'Cepheids.csv', sep=',')
        DF_dict['Cepheids'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I']
        DF_dict['Cepheids_anchors'] = pd.read_csv(data_tmp + 'Cepheids_anchors.csv', sep=',')
        DF_dict['Cepheids_anchors'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'mu', 'sig_mu']
        if fp.include_MW==True:
            print('Loading MW Cepheids...')
            DF_dict['Cepheids_MW'] = pd.read_csv(data_tmp+'Cepheids_MW.csv', sep=',')
            DF_dict['Cepheids_MW'].columns = ['Gal', 'logP', 'mW', 'sig_mW', 'M/H', 'z', 'V-I', 'pi', 'sig_pi']
        print('Loading SNe for the Cepheids-host galaxies...')
        DF_dict['SNe_Cepheids'] = pd.read_csv(data_tmp + 'SNe_Cepheids.csv', sep=',')
        DF_dict['SNe_Cepheids'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the TRGB related DataFrame
    if fp.include_TRGB==True:
        print('Loading TRGB & TRGB anchors...')
        DF_dict['TRGB'] = pd.read_csv(data_tmp+'TRGB.csv', sep=',')
        DF_dict['TRGB'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z']
        DF_dict['TRGB_anchors'] = pd.read_csv(data_tmp + 'TRGB_anchors.csv', sep=',')
        DF_dict['TRGB_anchors'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z', 'mu', 'sig_mu']
        print('Loading SNe for the TRGB-host galaxies...')
        DF_dict['SNe_TRGB'] = pd.read_csv(data_tmp + 'SNe_TRGB.csv', sep=',')
        DF_dict['SNe_TRGB'].columns = ['Gal', 'mB', 'sig_mB']

    # Loads the SNe for the aB fit
    if fp.fit_aB==True:
        print('Loading SNe for the aB fit...')
        DF_dict['SNe_Hubble'] = pd.read_csv(data_tmp + 'SNe_Hubble.csv', sep=',')
        DF_dict['SNe_Hubble'].columns = ['name', 'mB', 'sig_mB', 'z', 'sig_z']
        filter = ((DF_dict['SNe_Hubble']['z'] > fp.z_min) & (DF_dict['SNe_Hubble']['z'] < fp.z_max))
        SNe_other = DF_dict['SNe_Hubble'][~filter]
        DF_dict['SNe_Hubble'] = DF_dict['SNe_Hubble'][filter].reset_index(drop=True)
    else:
        SNe_other = None

    # Loads the table for the K_corrections
    if ((fp.Kcorr_Cep) == True or (fp.Kcorr_TRGB == True)):
        print('Loading the K-correction table...')
        table_Kcorr = pd.read_csv(data_static+'H1PStars/Anderson2022/tablea1.dat', delimiter=r"\s+")
        table_Kcorr = table_Kcorr.sort_values(by=['Teff', 'logg', '[Fe/H]', 'z', 'E(B-V)']).reset_index(drop=True)
    else:
        table_Kcorr = None

    return DF_dict, SNe_other, table_Kcorr