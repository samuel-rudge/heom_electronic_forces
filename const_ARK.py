# ----------------------------------------------------------------------------
#
#           CALCULATING THE ADAPTIVE RUNGE-KUTTA COEFFICIENTS 
#               FOR A TIME-INDEPENDENT MASTER EQUATION
#
# ----------------------------------------------------------------------------
# 
# Samuel Rudge 
#
# USAGE -       d,dhat = const_ARK.dvecs('Type')
#
# INPUTS - 
#       Type                        String with value 'Fehlberg' or 'Dormand Prince' 
#                                   for the two available methods
#
# OUTPUTS - 
#       d                           Vector of values for solution s.t. 
#                                   rho_{n+1} = rho_{n} + sum_{i} dt^{i}*d(i)*L*rho_{n}
#       
#       dhat                        Vector of values for error checking solution rhohat, calculated
#                                   in the same manner, with error related to abs(rho - rhohat)

import numpy as np

def dvecs(Type):

    if 'Fehlberg' in Type:

        Butcher_Fehlberg = np.array([
            [0,0,0,0,0,0],
            [1/4,0,0,0,0,0],
            [3/32,9/32,0,0,0,0],
            [1932/2197,-7200/2197,7296/2197,0,0,0],
            [439/216,-8,3680/513,-845/4104,0,0],
            [-8/27,2,-3544/2565,1859/4104,-11/40,0]],dtype=float)                                           # Define 
        b_Fehlberg = np.array([25/216,0,1408/2565,2197/4104,-1/5,0],dtype=float)
        bhat_Fehlberg = np.array([16/135,0,6656/12825,28561/56430,-9/50,2/55],dtype=float)

        d_Fehlberg = np.zeros((1,6),dtype=float)
        d_Fehlberg[0,0] = np.sum(b_Fehlberg)
        d_Fehlberg[0,1] = np.sum(b_Fehlberg*np.sum(Butcher_Fehlberg,axis=1))
        dhat_Fehlberg = np.zeros((1,6),dtype=float)
        dhat_Fehlberg[0,0] = np.sum(bhat_Fehlberg)
        dhat_Fehlberg[0,1] = np.sum(bhat_Fehlberg*np.sum(Butcher_Fehlberg,axis=1))


        for itrl in range(2,6):
            d_Fehlberg_cont = 0
            dhat_Fehlberg_cont = 0
            for itri in range(itrl,6):
                sum_Fehlberg = 0
                for itrj in range(itrl-1,itri):
                    prod_Fehlberg = Butcher_Fehlberg[itri,itrj]
                    for itrk in range(itrj-(itrl-2),itrj):
                        prod_Fehlberg = prod_Fehlberg*Butcher_Fehlberg[itrk+1,itrk]
                    sum_Fehlberg = sum_Fehlberg + prod_Fehlberg
                d_Fehlberg_cont = d_Fehlberg_cont + b_Fehlberg[itri]*sum_Fehlberg
                dhat_Fehlberg_cont = dhat_Fehlberg_cont + bhat_Fehlberg[itri]*sum_Fehlberg
            d_Fehlberg[0,itrl] = d_Fehlberg_cont
            dhat_Fehlberg[0,itrl] = dhat_Fehlberg_cont

        d = d_Fehlberg
        dhat = dhat_Fehlberg

    elif 'Dormand Prince' in Type:

        Butcher_DoPri = np.array([
            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
            [1/5,0.0,0.0,0.0,0.0,0.0,0.0],
            [3/40,9/40,0.0,0.0,0.0,0.0,0.0],
            [44/45,-56/15,32/9,0.0,0.0,0.0,0.0],
            [19372/6561,-25360/2187,64448/6561,-212/729,0.0,0.0,0.0],
            [9017/3168,-355/33,46732/5247,49/176,-5103/18656,0.0,0.0],
            [35/382,0.0,500/1113,125/192,-2187/6784,11/84,0.0]],dtype=float)
        b_DoPri = np.array([35/382,0.0,500/1113,125/192,-2187/6784,11/84,0.0],dtype=float)
        bhat_DoPri = np.array([5179/57600,0.0,7571/16695,393/640,-92097/339200,187/2100,1/40],dtype=float)

        d_DoPri = np.zeros((1,7),dtype=float)
        d_DoPri[0,0] = np.sum(b_DoPri)
        d_DoPri[0,1] = np.sum(b_DoPri*np.sum(Butcher_DoPri,axis=1))
        dhat_DoPri = np.zeros((1,7),dtype=float)
        dhat_DoPri[0,0] = np.sum(bhat_DoPri)
        dhat_DoPri[0,1] = np.sum(bhat_DoPri*np.sum(Butcher_DoPri,axis=1))

        for itrl in range(2,7):
            d_DoPri_cont = 0
            dhat_DoPri_cont = 0
            for itri in range(itrl,7):
                sum_DoPri = 0
                for itrj in range(itrl-1,itri):
                    prod_DoPri = Butcher_DoPri[itri,itrj]
                    for itrk in range(itrj-(itrl-2),itrj):
                        prod_DoPri = prod_DoPri*Butcher_DoPri[itrk+1,itrk]
                    sum_DoPri = sum_DoPri + prod_DoPri
                d_DoPri_cont = d_DoPri_cont + b_DoPri[itri]*sum_DoPri
                dhat_DoPri_cont = dhat_DoPri_cont + bhat_DoPri[itri]*sum_DoPri
            d_DoPri[0,itrl] = d_DoPri_cont
            dhat_DoPri[0,itrl] = dhat_DoPri_cont

        d = d_DoPri
        dhat = dhat_DoPri
    
    return d,dhat