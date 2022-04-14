'''
This module contains a function that fits the distance ladder to the different dataframes given.
'''
import numpy as np
import fit_parameters as fp


def fit_distance_ladder(DF_dict, cov_matrix = np.array([])):
    '''
    This function return the raw restult from the distance ladder fit, with the parameters from the fit_parameters.py.

    Parameters
    ----------
    DF_dict : dict of pandas DataFrame
        Dictionary that contains the different pandas DataFrame that have to be fitted.
    cov_matrix : numpy matrix
        Covariance matrix that will be used by the fit function. If empty, the function will build the covariance matrix
        from the DF_dict (By default cov_matrix is empty).

    Returns
    -------
    y : numpy array
        Signal vector y from the model equation y = L*q that has been fitted by the function.
    q_dict : dict of numpy array
        A dictionary that contains the parameters result and their uncertainty for the vector q from the model equation
        y = L*q that has been fitted by the function. A few other parameters (H0, chi2/2) are also in this dictionnary.
    L : numpy matrix
        Design matrix L from the model equation y = L*q that has been fitted by the function.
    Sigma2 : numpy matrix
        Covariance matrix resulting from the fit.
    '''


    ### Useful functions
    def logczexpansion(z):
        return np.log10(fp.c * z * (1 + 1 / 2 * (1 - fp.q0) * z - 1 / 6 * (1 - fp.q0 - 3 * fp.q0 ** 2 + fp.j0) * z ** 2))

    ### We first create an empty signal vector y, design matrix L and parameters vector q and diag (for Sigma)
    [y, L, diag] = 3 * [np.empty(0)]
    q_string = []

    ### For the Cepheids
    if fp.include_Cepheids == True:
        ### Load the DF
        # Use the Cepheids DF and check for its consistency with the Fit_parameters.py file
        Cepheids = DF_dict['Cepheids']
        galaxies_Cep = Cepheids['Gal'].drop_duplicates().reset_index(drop=True)  # List of Cepheids-host galaxy
        if len(galaxies_Cep) != fp.N_galaxies_Cep:
            print('ERROR: make sure that the value for N_galaxies_Cep in the Fit_parameters.py corresponds to the number of the galaxy in the Cepheids.csv file!')
            return
        # The anchors Cepheids
        Cepheids_anchors = DF_dict['Cepheids_anchors']
        anchors_Cep = Cepheids_anchors['Gal'].drop_duplicates().reset_index(drop=True)  # List of Cepheids-host galaxy
        if len(anchors_Cep) != fp.N_anchors_Cep:
            print('ERROR: make sure that the value for N_anchors_Cep in the Fit_parameters.py corresponds to the number of the galaxy in the Cepheids_anchors.csv file!')
            return
        # The MW Cepheids
        if fp.include_MW == True:
            Cepheids_MW = DF_dict['Cepheids_MW']
            if fp.multiple_zp == True:
                list_MW = Cepheids_MW['Gal'].drop_duplicates().reset_index(drop=True) # List of MW dataset
        #  The SNe from the Cepheid-host galaxy
        SNe_Cepheids = DF_dict['SNe_Cepheids']

        ### Complete y & diag for the covariance matrix
        # Cepheids
        y = np.append(y, Cepheids['mW'])
        err = (Cepheids['sig_mW'] + fp.added_scatter) ** 2
        diag = np.append(diag, err)
        # Anchors
        y = np.append(y, Cepheids_anchors['mW'] - Cepheids_anchors['mu'])
        err = (Cepheids_anchors['sig_mW'] + fp.added_scatter) ** 2 + Cepheids_anchors['sig_mu'] ** 2
        diag = np.append(diag, err)
        # External constraints (Delta mu from anchors = 0)
        y = np.append(y, np.zeros([1, fp.N_anchors_Cep]))
        index = Cepheids_anchors['Gal'].drop_duplicates().index  # find index of the 1st appearance of the anchor galaxy in the DF
        err = Cepheids_anchors.loc[index, 'sig_mu'].values ** 2
        diag = np.append(diag, err)
        # Milky Way Cepheids
        if fp.include_MW == True:
            if fp.fixed_zp == True:
                if fp.multiple_zp == True:
                    print('ERROR! you cannot use a fixed zp with multiple zp.')
                    return
                y = np.append(y, Cepheids_MW['mW'] \
                              - 10 + 5 * np.log10(Cepheids_MW['pi']) \
                              + 5 / np.log(10) * fp.zp / Cepheids_MW['pi'])
                err = (Cepheids_MW['sig_mW'] + fp.added_scatter) ** 2 \
                      + ((5 / np.log(10) / Cepheids_MW['pi'] - \
                          5 / np.log(10) * fp.zp / Cepheids_MW['pi'] ** 2) * Cepheids_MW['sig_pi']) ** 2 \
                      + (5 / np.log(10) / Cepheids_MW['pi'] * fp.sig_zp) ** 2
                diag = np.append(diag, err)
            else:
                y = np.append(y, Cepheids_MW['mW'] \
                              - 10 + 5 * np.log10(Cepheids_MW['pi']))
                err = (Cepheids_MW['sig_mW'] + fp.added_scatter) ** 2 \
                      + (5 / np.log(10) / Cepheids_MW['pi'] * Cepheids_MW['sig_pi']) ** 2
                diag = np.append(diag, err)
        # Put Zw [M/H] on the signal side (y) if fixed Zw
        if fp.fixed_Zw == True:
            minus_MH = Cepheids['M/H']
            minus_MH = np.append(minus_MH, Cepheids_anchors['M/H'])
            minus_MH = np.append(minus_MH, np.zeros([1, fp.N_anchors_Cep]))
            if fp.include_MW == True:
                minus_MH = np.append(minus_MH, Cepheids_MW['M/H'])
            y = y - fp.Zw * minus_MH
            diag = diag + (fp.sig_Zw * minus_MH) ** 2
        # SNe from Cepheid-host galaxies
        y = np.append(y, SNe_Cepheids['mB'])
        err = SNe_Cepheids['sig_mB'] ** 2
        diag = np.append(diag, err)

        ### Complete q_string for future dict
        for gal in galaxies_Cep:
            q_string.append(f'mu_{gal}')
        for gal in anchors_Cep:
            q_string.append(f'Dmu_{gal}')
        q_string.append('MW')
        q_string.append('MB')
        if fp.PLR_break == True:
            q_string.append('bs')
            q_string.append('bl')
        else:
            q_string.append('b')
        if fp.PLR_break2 == True:
            q_string.append('bL')
        if fp.fixed_Zw == False:
            q_string.append('Zw')
        if ((fp.include_MW == True) and (fp.fixed_zp == False)):
            if fp.multiple_zp == True:
                for i in range(len(list_MW)):
                    q_string.append(f'zp{i+1}')
            else:
                q_string.append('zp')

        ### fill L
        # First create the matrix
        L = np.zeros([len(y), len(q_string)])
        #  Cepheids
        index_offset = 0
        for i in range(len(Cepheids)):
            # mu
            gal = Cepheids.loc[i, 'Gal']
            index = q_string.index(f'mu_{gal}')
            L[i + index_offset, index] = 1
            # MW
            index = q_string.index('MW')
            L[i + index_offset, index] = 1
            #  logP
            if fp.PLR_break == False:
                index = q_string.index('b')
            else:
                threshold = np.log10(fp.break_P)
                if Cepheids.loc[i, 'logP'] < threshold:
                    index = q_string.index('bs')
                else:
                    index = q_string.index('bl')
            threshold2 = np.log10(fp.break_P2)
            if (fp.PLR_break2 == True) and (Cepheids.loc[i, 'logP'] >= threshold2):
                L[i + index_offset, index] = np.log10(fp.break_P2) - np.log10(fp.break_P)  # For the P1 -> P2 offset
                index = q_string.index('bL')
                L[i + index_offset, index] = Cepheids.loc[i, 'logP'] - np.log10(fp.break_P2)
            else:
                L[i + index_offset, index] = Cepheids.loc[i, 'logP']-np.log10(fp.break_P)
            #  Zw
            if fp.fixed_Zw == False:
                index = q_string.index('Zw')
                L[i + index_offset, index] = Cepheids.loc[i, 'M/H']
        #  Anchors Cepheids
        index_offset += i + 1
        for i in range(len(Cepheids_anchors)):
            # Dmu
            gal = Cepheids_anchors.loc[i, 'Gal']
            index = q_string.index(f'Dmu_{gal}')
            L[i + index_offset, index] = 1
            # MW
            index = q_string.index('MW')
            L[i + index_offset, index] = 1
            #  logP
            if fp.PLR_break == False:
                index = q_string.index('b')
            else:
                threshold = np.log10(fp.break_P)
                if Cepheids_anchors.loc[i, 'logP'] < threshold:
                    index = q_string.index('bs')
                else:
                    index = q_string.index('bl')
            L[i + index_offset, index] = Cepheids_anchors.loc[i, 'logP']-np.log10(fp.break_P)
            #  Zw
            if fp.fixed_Zw == False:
                index = q_string.index('Zw')
                L[i + index_offset, index] = Cepheids_anchors.loc[i, 'M/H']
        # External constriant (Delta mu)
        index_offset += i + 1
        for i in range(len(anchors_Cep)):
            gal = anchors_Cep[i]
            # Dmu
            index = q_string.index(f'Dmu_{gal}')
            L[i + index_offset, index] = 1
        # MW Cepheids
        if fp.include_MW == True:
            index_offset += i + 1
            for i in range(len(Cepheids_MW)):
                #  MW
                index = q_string.index('MW')
                L[i + index_offset, index] = 1
                #  logP
                if fp.PLR_break == False:
                    index = q_string.index('b')
                else:
                    threshold = np.log10(fp.break_P)
                    if Cepheids_MW.loc[i, 'logP'] < threshold:
                        index = q_string.index('bs')
                    else:
                        index = q_string.index('bl')
                L[i + index_offset, index] = Cepheids_MW.loc[i, 'logP']-np.log10(fp.break_P)
                #  Zw
                if fp.fixed_Zw == False:
                    index = q_string.index('Zw')
                    L[i + index_offset, index] = Cepheids_MW.loc[i, 'M/H']
                #  zp
                if (fp.fixed_zp == False):
                    if (fp.multiple_zp == False):
                        index = q_string.index('zp')
                        L[i + index_offset, index] = - 5 / np.log(10) / Cepheids_MW.loc[i, 'pi']
                    else:
                        index = list_MW[list_MW == Cepheids_MW.loc[i, 'Gal']].index[0]
                        index = q_string.index(f'zp{index+1}')
                        L[i + index_offset, index] = - 5 / np.log(10) / Cepheids_MW.loc[i, 'pi']

        #  SNe from Cepheid-host galaxies
        index_offset += i + 1
        for i in range(len(SNe_Cepheids)):
            # mu
            gal = SNe_Cepheids.loc[i, 'Gal']
            index = q_string.index(f'mu_{gal}')
            L[i + index_offset, index] = 1
            #  MB
            index = q_string.index('MB')
            L[i + index_offset, index] = 1

    ###  For the TRGB
    if fp.include_TRGB == True:
        ### Load the DF
        # TRGB
        TRGB = DF_dict['TRGB']
        galaxies_TRGB = TRGB['Gal'].drop_duplicates().reset_index(drop=True)  # List of TRGB-host galaxies
        if len(galaxies_TRGB) != fp.N_galaxies_TRGB:
            print('ERROR: make sure that the value for N_galaxies_TRGB in the Fit_parameters.py corresponds to \
                   the number of the galaxy in the TRGB.csv file!')
        #  TRGB anchors
        TRGB_anchors = DF_dict['TRGB_anchors']
        anchors_TRGB = TRGB_anchors['Gal'].drop_duplicates().reset_index(drop=True)  # List of TRGB-host galaxies
        if len(anchors_TRGB) != fp.N_anchors_TRGB:
            print('ERROR: make sure that the value for N_anchors_TRGB in the Fit_parameters.py corresponds to \
                   the number of the galaxy in the TRGB_anchors.csv file!')
        #  SNe_TRGB
        SNe_TRGB = DF_dict['SNe_TRGB']

        ### Complete y
        #  TRGB
        add_y = np.array(TRGB['m'] - TRGB['A'])
        err = TRGB['sig_m'] ** 2 + (0.5 * TRGB['A']) ** 2
        diag = np.append(diag, err)
        #  TRGB anchors
        add_y = np.append(add_y, TRGB_anchors['m'] - TRGB_anchors['A'] - TRGB_anchors['mu'])
        err = TRGB_anchors['sig_m'] ** 2 + (0.5 * TRGB_anchors['A']) ** 2 + TRGB_anchors['sig_mu'] ** 2
        diag = np.append(diag, err)
        if fp.use_color == True:
            color = TRGB['V-I'] - fp.mid_VI
            color = np.append(color, TRGB_anchors['V-I'] - fp.mid_VI)
            add_y = add_y - 0.2 * color
        y = np.append(y, add_y)
        # External constraint (Delta mu)
        y = np.append(y, np.zeros([1, fp.N_anchors_TRGB]))  # External constrainte (Delta mu_anchors = 0)
        index = TRGB_anchors[
            'Gal'].drop_duplicates().index  # find index of the 1st appearance of the anchor galaxy in the DF
        err = TRGB_anchors.loc[index, 'sig_mu'].values ** 2
        diag = np.append(diag, err)
        #  SNe TRGB
        y = np.append(y, SNe_TRGB['mB'])
        err = SNe_TRGB['sig_mB'] ** 2
        diag = np.append(diag, err)

        ### Complete q
        if fp.include_Cepheids == False:
            for gal in galaxies_TRGB:
                q_string.append(f'mu_{gal}')
            for gal in anchors_TRGB:
                q_string.append(f'Dmu_{gal}')
            q_string.append('MTRGB')
            q_string.append('MB')  #  Have to declare it here if there is no Cepheid
        else:
            for gal in galaxies_TRGB:
                if gal in list(galaxies_Cep):
                    if fp.different_mu == True:
                        q_string.append(f'Dmu_{gal}')  #  Delta mu if in both sets of galaxies
                else:
                    q_string.append(f'mu_{gal}')
            for gal in anchors_TRGB:
                if gal not in list(anchors_TRGB):
                    q_string.append(f'Dmu_{gal}')
            q_string.append('MTRGB')

        ### Complete L
        to_add = np.zeros([len(TRGB) + len(TRGB_anchors) + fp.N_anchors_TRGB + len(SNe_TRGB), len(q_string)])
        if fp.include_Cepheids == False:
            index_offset = 0
            # TRGB
            for i in range(len(TRGB)):
                # mu
                gal = TRGB.loc[i, 'Gal']
                index = q_string.index(f'mu_{gal}')
                to_add[i + index_offset, index] = 1
                # MTRGB
                index = q_string.index('MTRGB')
                to_add[i + index_offset, index] = 1
            index_offset += i + 1
            # TRGB anchors
            for i in range(len(TRGB_anchors)):
                # Dmu
                gal = TRGB_anchors.loc[i, 'Gal']
                index = q_string.index(f'Dmu_{gal}')
                to_add[i + index_offset, index] = 1
                # MTRGB
                index = q_string.index('MTRGB')
                to_add[i + index_offset, index] = 1
            index_offset += i + 1
            # External constriant (Delta mu)
            for i in range(len(TRGB_anchors)):
                # Dmu
                gal = TRGB_anchors.loc[i, 'Gal']
                index = q_string.index(f'Dmu_{gal}')
                to_add[i + index_offset, index] = 1
            index_offset += i + 1
            # SNe from TRGB-host galaxies
            for i in range(len(SNe_TRGB)):
                # mu
                index = galaxies_TRGB[galaxies_TRGB == SNe_TRGB.loc[i, 'Gal']].index[0]
                to_add[i + index_offset, index] = 1
                #  MB
                index = q_string.index('MB')
                to_add[i + index_offset, index] = 1
            L = to_add

        else:
            index_offset = 0
            # TRGB
            for i in range(len(TRGB)):
                # mu & Dmu
                gal = TRGB.loc[i, 'Gal']
                if gal in list(galaxies_TRGB) and gal in list(galaxies_Cep):
                    index = q_string.index(f'mu_{gal}')
                    to_add[i + index_offset, index] = 1
                    if fp.different_mu == True:
                        index = q_string.index(f'Dmu_{gal}')
                        to_add[i + index_offset, index] = 1
                else:
                    index = q_string.index(f'mu_{gal}')
                    to_add[i + index_offset, index] = 1
                    # MTRGB
                index = q_string.index('MTRGB')
                to_add[i + index_offset, index] = 1
            index_offset += i + 1
            # TRGB anchors
            for i in range(len(TRGB_anchors)):
                # Dmu
                gal = TRGB_anchors.loc[i, 'Gal']
                index = q_string.index(f'Dmu_{gal}')
                to_add[i + index_offset, index] = 1
                # MTRGB
                index = q_string.index('MTRGB')
                to_add[i + index_offset, index] = 1
            index_offset += i + 1
            # External constriant (Delta mu)
            for i in range(len(anchors_TRGB)):
                gal = anchors_TRGB[i]
                # Dmu
                index = q_string.index(f'Dmu_{gal}')
                to_add[i + index_offset, index] = 1
            # SNe from TRGB-host galaxies
            index_offset += i + 1
            for i in range(len(SNe_TRGB)):
                # mu & Dmu
                gal = SNe_TRGB.loc[i, 'Gal']
                if gal in list(galaxies_TRGB) and gal in list(galaxies_Cep):
                    index = q_string.index(f'mu_{gal}')
                    to_add[i + index_offset, index] = 1
                    if fp.different_mu == True:
                        index = q_string.index(f'Dmu_{gal}')
                        to_add[i + index_offset, index] = 1
                else:
                    index = q_string.index(f'mu_{gal}')
                    to_add[i + index_offset, index] = 1
                    #  MB
                index = q_string.index('MB')
                to_add[i + index_offset, index] = 1

            # Add a few lines to the right of the L matrix for the added parameters
            missing_columns = to_add.shape[1] - L.shape[1]
            L = np.hstack([L, np.zeros([len(L), missing_columns])])
            L = np.vstack([L, to_add])

    ###  For the aB fit
    if fp.fit_aB == True:
        ### Load the DF
        SNe_Hubble = DF_dict['SNe_Hubble']
        SNe_Hubble = SNe_Hubble[(SNe_Hubble['z'] > fp.z_min) & (SNe_Hubble['z'] < fp.z_max)].reset_index(drop=True)

        ### Complete y
        y = np.append(y, SNe_Hubble['mB'] - 5 * logczexpansion(SNe_Hubble['z']) - 25)
        err = SNe_Hubble['sig_mB'] ** 2
        diag = np.append(diag, err)

        ### Complete q
        q_string.append('5logH0')

        ### Complete L
        #  First add a column for the 1 extra parameters to the already existing matrice L
        L = np.hstack([L, np.zeros([len(L), 1])])
        #  Then add the lines
        to_add = np.zeros([len(SNe_Hubble), len(q_string)])
        #  MB columns
        to_add[:, q_string.index('MB')] = 1
        #  logH0 columns
        to_add[:, q_string.index('5logH0')] = -1
        L = np.vstack([L, to_add])

    ### Find the optimal parameters :
    if len(cov_matrix) == 0:
        C = np.diag(diag)
    else:
        C = cov_matrix
    LT = np.transpose(L)
    C1 = np.linalg.inv(C)
    Sigma2 = np.linalg.inv(np.matmul(np.matmul(LT, C1), L))
    q = np.matmul(np.matmul(np.matmul(Sigma2, LT), C1), y)
    chi2 = np.matmul(np.matmul(np.transpose(y - np.matmul(L, q)), C1), y - np.matmul(L, q))
    dof = len(y) - len(q)
    chi2_reduced = chi2 / dof

    # The q_dict now contains string name, q estimate, q uncertainty of the estimate
    sigma = np.sqrt(np.diag(Sigma2))
    q_dict = dict(zip(q_string, np.transpose([q,sigma])))

    # Complete the q_dict with H0 and the chi2_reduced
    ### Comlete the results
    # Compute the final value for H0 and add it to q_dict
    if fp.fit_aB == True:
        # finale value for H0 and
        [x, dx] = q_dict['5logH0']
        H0 = 10 ** (x / 5)
        sig_H0 = H0 * np.log(10) / 5 * dx
        q_dict['H0'] = [H0, sig_H0]
    else:
        # Compute H0 with the a_B from the Fit_parameters.py file
        [MB, sig_MB] = q_dict['MB']
        logH0 = 0.2 * MB + fp.aB + 5
        sig_logH0 = np.sqrt((0.2 * sig_MB) ** 2 + fp.sig_aB ** 2)
        H0 = 10 ** logH0
        sig_H0 = H0 * np.log(10) * sig_logH0
        q_dict['H0'] = [H0, sig_H0]
    # Add the chi2/dof
    q_dict['chi2/dof'] = [chi2_reduced, 0]

    return y, q_dict, L, Sigma2
