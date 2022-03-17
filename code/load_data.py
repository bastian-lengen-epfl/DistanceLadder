"""
This module contains a function that loads the previously pre-processed data and return the corresponding DataFrames
that are used by the other functions.
"""
import pandas as pd
import fit_parameters as fp
from astropy.io import fits
from astropy.table import Table
import warnings
from astropy.utils.exceptions import AstropyWarning

def load_data(data_tmp='./data_tmp/', data_static='./data_static/'):
    '''
    Loads and return DF_dict, a dictionary of pandas DataFrame of the different .csv. It also return the SNe that are
    out of the [z_min, z_max] interval in a SNe_other DataFrame

    :type   data_tmp: string
    :param  data_tmp: data_tmp directory, by default ./data_tmp/
    :type   data_static: string
    :param  data_static: data_static directory, by default ./data_static/
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

    # Remove warning about upper and lower-case management
    warnings.simplefilter('ignore', category=AstropyWarning)
    # Loads the table for the K_corrections
    if fp.Kcorr_Cep == True:
        print('Loading the Cepheids K-correction table...')
        with fits.open(data_static+'H1PStars/Anderson2022/tablea1.fits') as hdul:
            table_Kcorr_Cep = Table(hdul[1].data).to_pandas()
        table_Kcorr_Cep = table_Kcorr_Cep.sort_values(by=['Teff', 'logg', '[Fe/H]', 'z', 'E(B-V)'])\
                                         .reset_index(drop=True)
    else:
        table_Kcorr_Cep = None
    if fp.Kcorr_TRGB == True:
        print('Loading the TRGB K-correction table...')
        with fits.open(data_static+'H1PStars/Anderson2022/tablea3.fits') as hdul:
            table_Kcorr_TRGB = Table(hdul[1].data).to_pandas()
        table_Kcorr_TRGB = table_Kcorr_TRGB.sort_values(by=['Teff', 'logg', '[Fe/H]', 'z', 'E(B-V)'])\
                                           .reset_index(drop=True)
    else:
        table_Kcorr_TRGB = None

    return DF_dict, SNe_other, table_Kcorr_Cep, table_Kcorr_TRGB