################################################################
#                                                              #
#   COEFFICIENTS AND EXPONENTS OF BATH-CORRELATION FUNCTIONS   #
#                                                              # 
################################################################

import numpy as np
from matplotlib import pyplot as plt
# from input_parameters import *
# from constants import *
import barycentric_decomposition
import Pade_poles
import eta_gamma_pade

class bath_correlation_decomposition():

    def __init__(self,Ncutoff,specwidth,Nsupport_points_barycentric,Npoles_pade,symmetrized_fermi_specwidth,
                Temp,Nleads,Nsign,muvec,tol_Gamma_barycentric,tol_fermi_symmetrized_barycentric,
                wbl_YN,analytic_spectral_function_decomposition,tol_F):

        self.Ncutoff = Ncutoff
        self.specwidth = specwidth
        self.symmetrized_fermi_specwidth = symmetrized_fermi_specwidth
        self.Nsupport_points_barycentric = Nsupport_points_barycentric
        self.Npoles_pade = Npoles_pade
        self.Temp = Temp
        self.Nleads = Nleads
        self.Nsign = Nsign
        self.muvec = muvec
        self.tol_Gamma_barycentric = tol_Gamma_barycentric
        self.tol_fermi_symmetrized_barycentric = tol_fermi_symmetrized_barycentric
        self.wbl_YN = wbl_YN
        self.analytic_spectral_function_decomposition = analytic_spectral_function_decomposition
        self.tol_F = tol_F

        # ---------------------------------------------------------------
        #          AUTOMATICALLY GENERATE FREQUENCIES IN SPECTRUM
        # ---------------------------------------------------------------
        
        max_freq ="{:.2e}".format(self.Ncutoff*self.symmetrized_fermi_specwidth[0])                                 # Note: assuming same spectral density for each electrode - need to change if this ever changes in the future
        first_attempt = float(max_freq)*np.logspace(0,-10,num=int(self.Nsupport_points_barycentric/2))
        diff_vec = np.abs(first_attempt[1:] - first_attempt[:-1])
        greater_than_tol_F = len(diff_vec[diff_vec > self.tol_F])
        second_attempt = np.concatenate((first_attempt[0:greater_than_tol_F],np.flip(np.linspace(0,first_attempt[greater_than_tol_F],int(first_attempt[greater_than_tol_F]/self.tol_F)))),axis=0)
        self.Nsupport_points_barycentric = 2*len(second_attempt) - 1
        self.freq_vec_log = np.zeros((self.Nleads,self.Nsupport_points_barycentric),dtype=float)
        if self.wbl_YN == 1:
            for itrleads in range(self.Nleads):
                self.freq_vec_log[itrleads,:] = np.concatenate((-second_attempt,np.flip(second_attempt)[1:]),axis=0) + self.muvec[itrleads]
                # self.freq_vec_log[itrleads,:] = np.concatenate((-float(max_freq)*np.logspace(0,-3,num=int(self.Nsupport_points_barycentric/2))
                #                                  ,np.array([0],dtype=float),float(max_freq)*np.logspace(-3,0,num=int(self.Nsupport_points_barycentric/2))),axis=0) + self.muvec[itrleads]
        else:
            for itrleads in range(self.Nleads):
                self.freq_vec_log[itrleads,:] = np.concatenate((-second_attempt,np.flip(second_attempt)[1:]),axis=0) + self.muvec[itrleads]
                # self.freq_vec_log[itrleads,:] = np.concatenate((-float(max_freq)*np.logspace(0,-6,num=int(self.Nsupport_points_barycentric/2))
                #                                 ,np.array([0],dtype=float),float(max_freq)*np.logspace(-6,0,num=int(self.Nsupport_points_barycentric/2))),axis=0) + self.muvec[itrleads]

        self.spectral_functions()
        self.pade_decomposition()
        self.fermi_symmetrized_exact()
        self.barycentric_decomposition_gamma()
        self.barycentric_decomposition_fermi_symmetrized()

    # ---------------------------------------------------------------
    #                     SPECTRAL DECOMPOSITION 
    # ---------------------------------------------------------------

    def spectral_functions(self):

        self.Gamma_exact_log = np.zeros((self.Nleads,self.Nsupport_points_barycentric),dtype=float)
        if (self.wbl_YN == 0):
            for itrleads in range(self.Nleads):
                self.Gamma_exact_log[itrleads,:] = 2*np.pi*(self.specwidth[itrleads]**2)/(self.specwidth[itrleads]**2 \
                                                    + (self.freq_vec_log[itrleads,:] - self.muvec[itrleads])**2)

    def pade_decomposition(self):

        self.Pade = Pade_poles.Pade_spec_dec(self.Npoles_pade)
        self.Kappa,self.Zeta = self.Pade.get_Pade_parameters()

    def fermi_symmetrized_exact(self):
        
        self.fermi_symmetrized_exact_log = np.zeros((self.Nleads,self.Nsupport_points_barycentric),dtype=float)
        for itrleads in range(self.Nleads):
            self.fermi_symmetrized_exact_log[itrleads,:] = self.Pade.fermi_act((self.freq_vec_log[itrleads,:]-self.muvec[itrleads]),self.Temp) - 0.5

    def barycentric_decomposition_gamma(self):
        
        if self.analytic_spectral_function_decomposition == True:
            print("Spectral function (e.g. Lorentzian or in the WBL) can be analytically decomposed. Residues and Poles will be incorported \n \
                    directly into the gamma and eta matrices in the same way as for the Pade decomposition")
        else:
            mmax_while = 1
            all_imag_yn = True
            rel_err = 1.0
            real_poles_counter = 0
            while (rel_err > self.tol_Gamma_barycentric) or (all_imag_yn == True):
                mmax_while += 2
                BarycentricRational_Gamma,errors_barycentric = barycentric_decomposition.aaa(self.Gamma_exact_log[0,:],self.freq_vec_log[0,:],
                                                                                        mmax=mmax_while)
                rel_err_Gamma = errors_barycentric[-1]    
                poles_barycentric_temp,residues_barycentric_temp = BarycentricRational_Gamma.polres()
                npoles_decomposition = len(poles_barycentric_temp[np.imag(poles_barycentric_temp) > 0.0])
                all_imag_yn = any(np.imag(poles_barycentric_temp) == 0.0)
                if (all_imag_yn == True):
                    real_poles_counter += 1
                elif (all_imag_yn == False):
                    real_poles_counter = 0
                    largest_mmax_with_imag_component = mmax_while
                    smallest_rel_err_with_imag_component = rel_err_Gamma
                # print(all_imag_yn)
                # print(npoles_decomposition)
                # print(rel_err)
                if (real_poles_counter >= 5):
                    break
            print("Best barycentric decomposition of Gamma function is with "+str(int((largest_mmax_with_imag_component - 1)/2))+" barycentric poles and an error of " \
                    +str(smallest_rel_err_with_imag_component))
            self.Gamma_barycentric_log = np.zeros((self.Nleads,self.Nsupport_points_barycentric),dtype=float)
            self.Gamma_barycentric_err = np.zeros((self.Nleads),dtype=float)
            self.Npoles_Gamma_barycentric = np.zeros(self.Nleads,dtype=int)
            self.poles_Gamma_barycentric = np.zeros((self.Nleads,int((largest_mmax_with_imag_component - 1)/2),2),dtype=complex)
            self.residues_Gamma_barycentric = np.zeros((self.Nleads,int((largest_mmax_with_imag_component - 1)/2),2),dtype=complex)
            self.BR_Gamma = []
            for itrleads in range(self.Nleads):
                BarycentricRational_Gamma,errors_barycentric = barycentric_decomposition.aaa(self.Gamma_exact_log[itrleads,:],self.freq_vec_log[itrleads,:],
                                                                                            mmax=largest_mmax_with_imag_component)
                self.Gamma_barycentric_err[itrleads] = errors_barycentric[-1]
                poles_barycentric_temp,residues_barycentric_temp = BarycentricRational_Gamma.polres()
                self.Npoles_Gamma_barycentric[itrleads] = len(poles_barycentric_temp[np.imag(poles_barycentric_temp) > 0.0])
                for itrsign in range(2):
                    self.residues_Gamma_barycentric[itrleads,:,itrsign] = residues_barycentric_temp[(((-1)**(itrsign))*np.imag(poles_barycentric_temp)) > 0.0]
                    self.poles_Gamma_barycentric[itrleads,:,itrsign] = poles_barycentric_temp[(((-1)**(itrsign))*np.imag(poles_barycentric_temp)) > 0.0]
                    self.residues_Gamma_barycentric[itrleads,:,itrsign] = \
                        self.residues_Gamma_barycentric[itrleads,np.abs(np.imag(self.poles_Gamma_barycentric[itrleads,:,itrsign])).argsort(),itrsign]
                    self.poles_Gamma_barycentric[itrleads,:,itrsign] = \
                        self.poles_Gamma_barycentric[itrleads,np.abs(np.imag(self.poles_Gamma_barycentric[itrleads,:,itrsign])).argsort(),itrsign]
                self.Gamma_barycentric_log[itrleads,:] = BarycentricRational_Gamma(self.freq_vec_log[itrleads,:])
                self.BR_Gamma.append(BarycentricRational_Gamma)

    def barycentric_decomposition_fermi_symmetrized(self):
        
        mmax_while = 1
        all_imag_yn = True
        rel_err_fermi_symmetrized = 1.0
        real_poles_counter = 0
        while (rel_err_fermi_symmetrized > self.tol_fermi_symmetrized_barycentric) or (all_imag_yn == True):
            mmax_while += 2
            BarycentricRational_fermi_symmetrized,errors_barycentric = barycentric_decomposition.aaa(self.fermi_symmetrized_exact_log[0,:],self.freq_vec_log[0,:],
                                                                                    mmax=mmax_while)
            rel_err_fermi_symmetrized = errors_barycentric[-1]    
            poles_barycentric_temp,residues_barycentric_temp = BarycentricRational_fermi_symmetrized.polres()
            npoles_decomposition = len(poles_barycentric_temp[np.imag(poles_barycentric_temp) > 0.0])
            all_imag_yn = any(np.imag(poles_barycentric_temp) == 0.0)
            if (all_imag_yn == True):
                real_poles_counter += 1
            elif (all_imag_yn == False):
                real_poles_counter = 0
                largest_mmax_with_imag_component = mmax_while
                smallest_rel_err_with_imag_component = rel_err_fermi_symmetrized
            # print(all_imag_yn)
            # print(npoles_decomposition)
            # print(rel_err_fermi_symmetrized)
            if (real_poles_counter >= 10):
                break
        print("Best barycentric decomposition of symmetrized Fermi-Dirac function is with "+str(int((largest_mmax_with_imag_component - 1)/2))+" barycentric poles and an error of " \
                +str(smallest_rel_err_with_imag_component))
        self.fermi_symmetrized_barycentric_log = np.zeros((self.Nleads,self.Nsupport_points_barycentric),dtype=float)
        self.fermi_symmetrized_barycentric_err = np.zeros((self.Nleads),dtype=float)
        self.Npoles_fermi_symmetrized_barycentric = int((largest_mmax_with_imag_component - 1)/2)
        self.poles_fermi_symmetrized_barycentric = np.zeros((self.Nleads,int((largest_mmax_with_imag_component - 1)/2),2),dtype=complex)
        self.residues_fermi_symmetrized_barycentric = np.zeros((self.Nleads,int((largest_mmax_with_imag_component - 1)/2),2),dtype=complex)
        self.BR_fermi_symmetrized = []
        for itrleads in range(self.Nleads):
            BarycentricRational_fermi_symmetrized,errors_barycentric = barycentric_decomposition.aaa(self.fermi_symmetrized_exact_log[itrleads,:],self.freq_vec_log[itrleads,:],
                                                                                        mmax=largest_mmax_with_imag_component)
            self.fermi_symmetrized_barycentric_err[itrleads] = errors_barycentric[-1]
            poles_barycentric_temp,residues_barycentric_temp = BarycentricRational_fermi_symmetrized.polres()
            for itrsign in range(2):
                self.residues_fermi_symmetrized_barycentric[itrleads,:,itrsign] = residues_barycentric_temp[(((-1)**(itrsign))*np.imag(poles_barycentric_temp)) > 0.0]
                self.poles_fermi_symmetrized_barycentric[itrleads,:,itrsign] = poles_barycentric_temp[(((-1)**(itrsign))*np.imag(poles_barycentric_temp)) > 0.0]
                self.residues_fermi_symmetrized_barycentric[itrleads,:,itrsign] = \
                    self.residues_fermi_symmetrized_barycentric[itrleads,np.abs(np.imag(self.poles_fermi_symmetrized_barycentric[itrleads,:,itrsign])).argsort(),itrsign]
                self.poles_fermi_symmetrized_barycentric[itrleads,:,itrsign] = \
                    self.poles_fermi_symmetrized_barycentric[itrleads,np.abs(np.imag(self.poles_fermi_symmetrized_barycentric[itrleads,:,itrsign])).argsort(),itrsign]
            self.fermi_symmetrized_barycentric_log[itrleads,:] = BarycentricRational_fermi_symmetrized(self.freq_vec_log[itrleads,:])
            self.BR_fermi_symmetrized.append(BarycentricRational_fermi_symmetrized)

    # ---------------------------------------------------------------
    #                   BATH-CORRELATION SPECTRUMS
    # ---------------------------------------------------------------

    def return_approximated_functions_exact(self):

        return self.Gamma_exact_log,self.fermi_symmetrized_exact_log

    def return_approximated_functions_barycentric(self):

        if self.analytic_spectral_function_decomposition == False:
            return self.Gamma_barycentric_log,self.fermi_symmetrized_barycentric_log
        else:
            return self.fermi_symmetrized_barycentric_log

    # ---------------------------------------------------------------
    #             BARYCENTRIC BATH-CORRELATION EXPANSION
    # ---------------------------------------------------------------

    def barycentric_bath_correlation_expansion(self):

        if self.analytic_spectral_function_decomposition == True:
            if self.wbl_YN == 0:
                self.eta_vec_barycentric = np.zeros((self.Nleads,self.Nsign,self.Npoles_fermi_symmetrized_barycentric+1),dtype=complex)
                self.gamma_vec_barycentric = np.zeros((self.Nleads,self.Nsign,self.Npoles_fermi_symmetrized_barycentric+1),dtype=complex)
                for itrleads in range(self.Nleads):
                    for itrsign in range(2):
                        self.eta_vec_barycentric[itrleads,itrsign,0] = np.pi*self.specwidth[itrleads]*(0.5 + ((-1)**itrsign)*
                                                                        self.BR_fermi_symmetrized[itrleads](self.muvec[itrleads] + ((-1)**itrsign)*1j*self.specwidth[itrleads]))
                        self.gamma_vec_barycentric[itrleads,itrsign,0] = self.specwidth[itrleads] - ((-1)**itrsign)*1j*self.muvec[itrleads]
                        for itrpole in range(self.Npoles_fermi_symmetrized_barycentric):
                            self.gamma_vec_barycentric[itrleads,itrsign,itrpole+1] = -1j*((-1)**itrsign)*self.poles_fermi_symmetrized_barycentric[itrleads,itrpole,itrsign]
                            self.eta_vec_barycentric[itrleads,itrsign,itrpole+1] = 1j*self.residues_fermi_symmetrized_barycentric[itrleads,itrpole,itrsign]*\
                                                        self.spectral_function_Lorentzian(self.poles_fermi_symmetrized_barycentric[itrleads,itrpole,itrsign],\
                                                            self.specwidth[itrleads],self.muvec[itrleads])
            else:
                self.eta_vec_barycentric = np.zeros((self.Nleads,self.Nsign,self.Npoles_fermi_symmetrized_barycentric),dtype=complex)
                self.gamma_vec_barycentric = np.zeros((self.Nleads,self.Nsign,self.Npoles_fermi_symmetrized_barycentric),dtype=complex)
                for itrleads in range(self.Nleads):
                    for itrsign in range(2):
                        for itrpole in range(self.Npoles_fermi_symmetrized_barycentric):
                            self.gamma_vec_barycentric[itrleads,itrsign,itrpole] = -1j*((-1)**itrsign)*self.poles_fermi_symmetrized_barycentric[itrleads,itrpole,itrsign]
                            self.eta_vec_barycentric[itrleads,itrsign,itrpole] = 1j*self.residues_fermi_symmetrized_barycentric[itrleads,itrpole,itrsign]*2*np.pi

        return self.eta_vec_barycentric,self.gamma_vec_barycentric
        # return self.poles_barycentric_temp,self.residues_barycentric_temp,self.real_poles,self.barycentric_err

    # ---------------------------------------------------------------
    #                 PADE BATH-CORRELATION EXPANSION
    # ---------------------------------------------------------------

    def pade_bath_correlation_expansion(self):

        if self.wbl_YN == 0:
            self.eta_vec_pade,self.gamma_vec_pade = eta_gamma_pade.genbathcoeff(kappa=self.Kappa,zeta=self.Zeta,specwidth=self.specwidth,muvec=self.muvec,
                                                temp=self.Temp,nleads=self.Nleads,npoles=self.Npoles_pade,nsign=self.Nsign)
        elif self.wbl_YN == 1:
            self.eta_vec_pade,self.gamma_vec_pade = eta_gamma_pade.genbathcoeff_wbl(kappa=self.Kappa,zeta=self.Zeta,muvec=self.muvec,
                                                temp=self.Temp,nleads=self.Nleads,npoles=self.Npoles_pade,nsign=self.Nsign)

        return self.eta_vec_pade,self.gamma_vec_pade
                
    def return_freq_vec(self):
        return self.freq_vec_log,self.Nsupport_points_barycentric

    def spectral_function_Lorentzian(self,w,specwidth,mu):

        g = 2*np.pi*(specwidth**2)/(specwidth**2 + (w - mu)**2)
        return g

    def fermi_symmetrized(self,w,Temp,mu):

        f = 1/(np.exp((w-mu)/Temp)+1) - 0.5
        return f

