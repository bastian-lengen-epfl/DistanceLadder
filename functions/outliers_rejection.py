"""
This module contains the outliers rejections process that will be called by the pipeline.
"""
import numpy as np
import Fit_parameters as Fp
from Fitting import fit_distance_ladder


def single_kappa_clipping(DF_dict, kappa=2.7):
    '''
    Return both the DF_dict after the outlier rejection and the DF_dict which will contain all the outliers
    rejected from the DF_dict by following the single kappa clipping algorithm. This algorithm consists of
    rejecting the worst point, one at a time. The fit is done. If the worst point is at more than kappa*std(errors)
    from the fit, it is rejected. After 1 point is rejected, the procedure restart. The fit is done and the
    worst point is rejected at each iteration. The process stop when no point is worse than kappa*std(errors).
    Note that the algorithm only reject points from the Cepheids, Cepheids_anchors, Cepheids_MW or SNe Hubble

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that will be fitted.
    :type   kappa: float
    :param  kappa: The kappa parameters for which the algorithm has to stop.
    '''

    ### Create the DF_dict_outliers for the outliers and load the DF_dict
    DF_dict_outliers = dict()
    if Fp.include_Cepheids == True:
        DF_dict_outliers['Cepheids'] = []
        DF_dict_outliers['Cepheids_anchors'] = []
        Cepheids = DF_dict['Cepheids']
        Cepheids_anchors = DF_dict['Cepheids_anchors']
        N_Cep = len(Cepheids)
        N_anc = len(Cepheids_anchors)
        N_MW = 0  # Set to 0 if there is no MW cepheids
        if Fp.include_MW == True:
            DF_dict_outliers['Cepheids_MW'] = []
            Cepheids_MW = DF_dict['Cepheids_MW']
            N_MW = len(Cepheids_MW)
    if Fp.fit_aB == True:
        DF_dict_outliers['SNe_Hubble'] = []
        SNe_Hubble = DF_dict['SNe_Hubble']
        N_SN = len(SNe_Hubble)


    ### First iteration
    y, q_dict, L = fit_distance_ladder(DF_dict)
    q = np.array([]) # Load the q values
    for str in q_dict:
        q = np.append(q_dict[str][0])

    ### Initialize the algorithm
    errors = np.abs(np.array(y - np.matmul(L, q)))
    std = np.std(errors)
    worst_Cep, worst_SNe = 0, 0 # Set to 0 if there is not both of them for the later comparison
    if Fp.include_Cepheids == True:
        worst_Cep = np.max(errors[:N_Cep+N_anc+N_MW])
    if Fp.fit_aB == True:
        worst_SNe = np.max(errors[-N_SN:])
    # Compare them
    worst = np.max(worst_SNe, worst_Cep)

    ### Start iterations
    while worst >= kappa*std:
        index_worst = list(errors).index(worst)
        # Get the outlier in the DF_dict_outliers and delete it from the DF_dict
        if Fp.include_Cepheids == True:
            if index_worst<N_Cep:
                index = index_worst
                DF_dict_outliers['Cepheids'] = DF_dict_outliers['Cepheids']\
                                             .append(Cepheids.iloc[index])
                Cepheids.drop(index=index).reset_index(drop=True)
            elif index_worst<N_Cep+N_anc:
                index = index_worst-N_Cep
                DF_dict_outliers['Cepheids_anchors'] = DF_dict_outliers['Cepheids_anchors']\
                                                    .append(Cepheids_anchors.iloc[index])
                Cepheids_anchors.drop(index=index).reset_index(drop=True)
            elif ((Fp.include_MW == True) and (index_worst<N_Cep+N_anc+N_MW)):
                index = index_worst-N_Cep-N_anc
                DF_dict_outliers['Cepheids_MW'] = DF_dict_outliers['Cepheids_MW']\
                                                .append(Cepheids_MW.iloc[index])
                Cepheids_MW.drop(index=index).reset_index(drop=True)
        if ((Fp.fit_aB == True) and (index_worst>N_Cep+N_anc+N_MW)):
            offset = len(y)-N_SN
            index = index + offset
            DF_dict_outliers['SNe_Hubble'] = DF_dict_outliers['SNe_Hubble']\
                                           .append(SNe_Hubble.iloc[index])
            SNe_Hubble.drop(index=index).reset_index(drop=True)

        # Re-iterate
        y, q_dict, L = fit_distance_ladder(DF_dict)
        q = np.array([])  # Load the q values
        for str in q_dict:
            q = np.append(q_dict[str][0])

        errors = np.abs(np.array(y - np.matmul(L, q)))
        std = np.std(errors)
        worst_Cep, worst_SNe = 0, 0  # Set to 0 if there is not both of them for the later comparison
        if Fp.include_Cepheids == True:
            worst_Cep = np.max(errors[:N_Cep + N_anc + N_MW])
        if Fp.fit_aB == True:
            worst_SNe = np.max(errors[-N_SN:])
        # Compare them
        worst = np.max(worst_SNe, worst_Cep)

    return DF_dict, DF_dict_outliers
