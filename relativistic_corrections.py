"""
This module contains the functions that correct for different relativistic effect by following the methodology
from Anderson (2019, 2022), [2019A&A...631A.165A] and [2022A&A...658A.148A].
"""
import numpy as np
import fit_parameters as fp

def RLB_correction(DF_dict):
    '''
    Return the corrected the DF_dict for the Redshift Leavitt Bias (RLB) following the methodology from
    Anderson (2019) [2019A&A...631A.165A] for Cepheids

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that will be fitted.
    '''
    Cepheids = DF_dict['Cepheids']
    Cepheids_anchors = DF_dict['Cepheids_anchors']

    Cepheids['logP'] = Cepheids['logP'] - np.log10(1+Cepheids['z'])
    Cepheids_anchors['logP'] = Cepheids_anchors['logP'] - np.log10(1 + Cepheids_anchors['z'])

    return DF_dict

def K_corr_Cep(DF_dict, filter='W'):
    '''
    Return the corrected the DF_dict for the Cepheids K-corrections following the methodology from Anderson (2021)
    [2022A&A...658A.148A]

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that will be fitted.
    :type   filter: string
    :param  filter: filter in which the K-corrections have to be applied, by default the wenseheit W
    '''

    # Defines an interpolate function to interpolate between reference points from Anderson (2021)
    def interpolate(z, z_ref, m_ref, c_ref):
        '''
        Returns the m and c of the linear interpolation between the reference values from Anderson (2022)

        :param z:       redshif of the Cepheid
        :type:          double
        :param z_ref:   Array containing all reference redshift from Anderson (2021)
        :type:          numpy array
        :param m_ref:   Array containing all reference slope m from Anderson (2021)
        :type:          numpy array
        :param c_ref:   Array containing all reference intercept c from Anderson (2021)
        :type:          numpy array
        '''
        # Linear interpolation
        for i in range(len(z_ref) - 1):
            if z_ref[i] <= z and z < z_ref[i + 1]:
                m = m_ref[i] + (z - z_ref[i]) * (m_ref[i + 1] - m_ref[i]) / (z_ref[i + 1] - z_ref[i])
                c = c_ref[i] + (z - z_ref[i]) * (c_ref[i + 1] - c_ref[i]) / (z_ref[i + 1] - z_ref[i])
                return m, c
            else:
                pass
        if z < z_ref[0]:
            m = m_ref[0]
            c = c_ref[0]
        else:
            m = m_ref[-2]
            c = c_ref[-2]

        return m, c

    # Reference points from Anderson (2021)
    z_ref = np.array([0.0019, 0.0056, 0.0098, 0.0172, 0.0245])
    if filter == 'W':
        m_ref = np.array([3.48, 2.68, 1.89, 1.07, 0.31]) * 1e-3
        c_ref = np.array([0.51, 1.74, 3.25, 5.96, 8.05]) * 1e-3
    elif filter == 'F555W':
        m_ref = np.array([-2.84, -8.65, -15.16, -26.85, -38.66]) * 1e-3
        c_ref = np.array([-1.74, -5.47, -9.48, -15.67, -20.51]) * 1e-3
    elif filter == 'F814W':
        m_ref = np.array([-1.02, -3.11, -5.47, -9.40, -12.73]) * 1e-3
        c_ref = np.array([-0.17, -0.91, -1.79, -2.82, -4.02]) * 1e-3
    elif filter == 'F160W':
        m_ref = np.array([-1.18, -3.53, -6.04, -10.10, -14.38]) * 1e-3
        c_ref = np.array([1.00, 1.93, 3.19, 5.69, 8.28]) * 1e-3
    else:
        print('Error, choose a filter in which the Cepheids K-corrections is available !')

    # Load the different DataFrames
    Cepheids = DF_dict['Cepheids']
    Cepheids_anchors = DF_dict['Cepheids_anchors']

    # Correct each DataFrame for it
    for i in Cepheids.index:
        m, c = interpolate(Cepheids.loc[i, 'z'], z_ref, m_ref, c_ref)
        Cepheids.loc[i, 'mW'] = Cepheids.loc[i, 'mW'] \
                              + (m * Cepheids.loc[i, 'logP'] + c) * Cepheids.loc[i,'z'] * Cepheids.loc[i,'V-I'] \
                              - 0.105 * Cepheids.loc[i,'z'] * Cepheids.loc[i,'V-I'] # For F99 redshift law
    for i in Cepheids_anchors.index:
        m, c = interpolate(Cepheids_anchors.loc[i, 'z'], z_ref, m_ref, c_ref)
        Cepheids_anchors.loc[i, 'mW'] = Cepheids_anchors.loc[i, 'mW'] \
                              + (m * Cepheids_anchors.loc[i, 'logP'] + c) * Cepheids_anchors.loc[i,'z'] * Cepheids_anchors.loc[i,'V-I'] \
                              - 0.105 * Cepheids_anchors.loc[i,'z'] * Cepheids_anchors.loc[i,'V-I'] # For F99 redshift law

    if fp.include_MW == True:
        Cepheids_MW = DF_dict['Cepheids_MW']
        for i in Cepheids_MW.index:
            m, c = interpolate(Cepheids_MW.loc[i, 'z'], z_ref, m_ref, c_ref)
            Cepheids_MW.loc[i, 'MW'] = Cepheids_MW.loc[i, 'mW'] \
                              + (m * Cepheids_MW.loc[i, 'logP'] + c) * Cepheids_MW.loc[i,'z'] * Cepheids_MW.loc[i,'V-I'] \
                              - 0.105 * Cepheids_MW.loc[i,'z'] * Cepheids_MW.loc[i,'V-I'] # For F99 redshift law

    return DF_dict

def K_corr_TRGB(DF_dict, filter='I'):
    '''
    Return the corrected the DF_dict for the TRGB K-corrections following the methodology from Anderson (2022)
    [2022A&A...658A.148A]

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that will be fitted.
    :type   filter: string
    :param  filter: filter in which the K-corrections have to be applied, by default the I-band (F814W)
    '''
    # Parameters from Anderson (2021)
    if filter == 'V':
        a, b = -0.0012, -4.1162
    elif filter == 'I':
        a, b = -0.0004, -1.4075
    elif filter == 'H':
        a, b = 0.0001, -1.6241
    else:
        print('Error, choose a filter in which the Cepheids K-corrections is available !')
        return

    # Load the different DataFrames
    TRGB = DF_dict['TRGB']
    TRGB_anchors = DF_dict['TRGB_anchors']

    # Apply the correction:
    for i in TRGB.index:
        TRGB.loc[i,'m'] = TRGB.loc[i,'m'] \
                        + (a + b * TRGB.loc[i,'z']) * TRGB.loc[i,'z'] * TRGB.loc[i,'V-I']
    for i in TRGB_anchors.index:
        TRGB_anchors.loc[i,'m'] = TRGB_anchors.loc[i,'m'] \
                                + (a + b * TRGB_anchors.loc[i,'z']) * TRGB_anchors.loc[i,'z'] * TRGB_anchors.loc[i,'V-I']
    return DF_dict
