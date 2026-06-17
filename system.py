import numpy as np
from numpy.core.numeric import ones
from constants import *  # pylint: disable=unused-import
from input_parameters import * #  pylint: disable=unused-import
import CreAnn
import Franck_Condon
import itertools

def system_operators(Single_El_Int,Double_El_Int,Vib_Freq_qu,El_Nuclear_Couplings_qu,
                     El_Nuclear_Couplings_cl,Nel,N_qu_vib_modes,N_cl_vib_modes,max_occ_qu_vib_modes,dim_rho,
                     len_x_vec,x_vec,dx_heom,small_polaron_yn,nondiag_key):

    """
    Generate all system related operators and associated Fock states
    """

    # Generate raw creation and annihilation operators from CreAnn class
    """
    For explanation of parameters, see input_parameters.py. For explanation of creation and annihilation operators, see CreAnn.py
    """

    Constraints = [Nel,N_qu_vib_modes,max_occ_qu_vib_modes] # Define constraints of system to be inputted into CreAnn function
    if bool(N_qu_vib_modes):
        CreAnn1 = CreAnn.CreAnn(Constraints,'Both') # Run program to generate creation and annihilation operators and Fock space
        d_ops,d,ddag,b_ops,b,bdag,Fock_states = CreAnn1.return_operators() # pylint: disable=unbalanced-tuple-unpacking
    else:
        CreAnn1 = CreAnn.CreAnn(Constraints,'Fermi') # Run program to generate creation and annihilation operators and Fock space
        d_ops,d,ddag,Fock_states = CreAnn1.return_operators() # pylint: disable=unbalanced-tuple-unpacking

    # Generate electronic part of Hamiltonian

    Ham_el = np.zeros((dim_rho,dim_rho),dtype=complex) 
    if Nel == 1:
        if bool(small_polaron_yn):
            if not nondiag_key:
                Ham_el += (Single_El_Int[0,0] - small_polaron_shift_qu)*np.matmul(ddag[:,:,0],d[:,:,0])
            else:
                Ham_el += (Single_El_Int[0,0])*np.matmul(ddag[:,:,0],d[:,:,0])    
        else:
            Ham_el += (Single_El_Int[0,0])*np.matmul(ddag[:,:,0],d[:,:,0])
    elif Nel > 1:
        for itr_el1 in range(Nel):
            for itr_el2 in range(Nel):
                Ham_el += (Single_El_Int[itr_el1,itr_el2])*np.matmul(ddag[:,:,itr_el1],d[:,:,itr_el2])
                if (itr_el1 < itr_el2):
                    Ham_el += Double_El_Int[itr_el1,itr_el2]*np.matmul(ddag[:,:,itr_el1],np.matmul(d[:,:,itr_el1],np.matmul(ddag[:,:,itr_el2],d[:,:,itr_el2])))

    Ham = Ham_el

    if N_qu_vib_modes > 0:

        # Generate quantum nuclear part of Hamiltonian
        
        Ham_vib_qu = np.zeros((dim_rho,dim_rho),dtype=complex)
        nvib_qu = np.zeros((dim_rho,dim_rho,N_qu_vib_modes),dtype=complex)
        for itr_qu_vib_modes in range(N_qu_vib_modes):
            nvib_qu[:,:,itr_qu_vib_modes] = np.matmul(bdag[:,:,itr_qu_vib_modes],b[:,:,itr_qu_vib_modes])
            Ham_vib_qu += Vib_Freq_qu[itr_qu_vib_modes]*nvib_qu[:,:,itr_qu_vib_modes]
        nvib_qu_diag = np.copy(nvib_qu)
        nvib_qu_nondiag = np.zeros((dim_rho,dim_rho,N_qu_vib_modes),dtype=float)
        nvib_qu_nondiag[:,:,0] = np.matmul(bdag[:,:,0] + El_Nuclear_Couplings_qu[0]/Vib_Freq_qu[0]*\
                                np.matmul(ddag[:,:,0],d[:,:,0]),\
                                b[:,:,0] + El_Nuclear_Couplings_qu[0]/Vib_Freq_qu[0]*\
                                np.matmul(ddag[:,:,0],d[:,:,0]))
        
        Ham += Ham_vib_qu

        # Generate Franck-Condon matrix and dressed fermionic annihilation and creation operators
        FC_object = Franck_Condon.Franck_Condon(Constraints,El_Nuclear_Couplings_qu,Vib_Freq_qu[0])
        FC_Matrix,FC_Matrix_Fock_Space = FC_object.return_FC_Operators()
        d_ops_dressed = np.zeros((dim_rho,dim_rho,Nel,2),dtype=float)
        for itr_el in range(Nel):
            d_ops_dressed[:,:,itr_el,0] = np.matmul(d_ops[:,:,itr_el,0],scipy.linalg.expm(El_Nuclear_Couplings_qu[0]/Vib_Freq_qu[0]*(bdag[:,:,0] - b[:,:,0])))
            d_ops_dressed[:,:,itr_el,1] = np.matmul(d_ops[:,:,itr_el,1],scipy.linalg.expm(-El_Nuclear_Couplings_qu[0]/Vib_Freq_qu[0]*(bdag[:,:,0] - b[:,:,0])))
        # Kmatrix = np.matmul(d[:,:,0],ddag[:,:,0])
        # ind_occ_states = np.where(Fock_states[:,0] == 1)
        # pairs_occ_states = list(zip(*list(itertools.product(ind_occ_states[0], repeat=2))))
        # Kmatrix[pairs_occ_states[0],pairs_occ_states[1]] = FC_Matrix[:,:,0].reshape(dim_vib_mode_qu**2)
        Kmatrix = scipy.linalg.expm(np.matmul(np.matmul(ddag[:,:,0],d[:,:,0]),
                                    El_Nuclear_Couplings_qu[0]/Vib_Freq_qu[0]*
                                    (bdag[:,:,0] - b[:,:,0])))
        # Kmatrix_inv = Kmatrix.transpose()
        Kmatrix_inv = scipy.linalg.inv(Kmatrix)
        
        if nondiag_key or not small_polaron_yn:
            # Generate electronic-vibrational interaction Hamiltonian
            Ham_qu_elvib_int = np.zeros((dim_rho,dim_rho),dtype=complex)
            for itr_qu_elvib_int in range(N_el_vib_int_qu):
                Ham_qu_elvib_int += El_Nuclear_Couplings_qu[itr_qu_elvib_int]*\
                                        np.matmul(np.matmul(ddag[:,:,0],d[:,:,0]),
                                        (bdag[:,:,itr_qu_elvib_int]+b[:,:,itr_qu_elvib_int]))
            Ham += Ham_qu_elvib_int
                        
    # Generate partial derivative of Hamiltonian with respect to semi-classical mode

    deriv_Ham = np.zeros((dim_rho,dim_rho,N_el_vib_int_cl,len_x_vec),dtype=complex)
    for itr_nuc_cl in range(N_cl_vib_modes):
        # for itr_el_vib_int_cl in range(N_el_vib_int_cl):
        for itrx in range(len_x_vec):
            deriv_Ham[:,:,itr_nuc_cl,itrx] = El_Nuclear_Couplings_cl[itr_nuc_cl]*(np.matmul(ddag[:,:,0],d[:,:,1])+
                                                                                  np.matmul(ddag[:,:,1],d[:,:,0]))

    # Generate q-dependent Hamiltonian

    Ham_x = np.zeros((dim_rho,dim_rho,len_x_vec,5),dtype=complex)
    dx_vec = np.array([0,2*dx_heom,dx_heom,-dx_heom,-2*dx_heom],dtype=complex)
    for itrx in range(len_x_vec):
        for itrdx in range(5):
            Ham_x[:,:,itrx,itrdx] = Ham + El_Nuclear_Couplings_cl[0]*(x_vec[itrx]+dx_vec[itrdx])*\
                                    (np.matmul(ddag[:,:,0],d[:,:,1])+np.matmul(ddag[:,:,1],d[:,:,0]))

    if not bool(N_qu_vib_modes):
        return d_ops,d,ddag,Fock_states,Ham,Ham_x,deriv_Ham
    else:
        return d_ops,d,ddag,Fock_states,Ham,Ham_x,deriv_Ham,b_ops,b,bdag,nvib_qu,nvib_qu_diag,nvib_qu_nondiag,\
                d_ops_dressed,FC_Matrix,FC_Matrix_Fock_Space,Kmatrix,Kmatrix_inv

if __name__=='__main__': # Example parameters if code is run directly from terminal and not called as a function.

    Nel = 2                                                                       
    N_qu_vib_modes = 1 # Must be 1 for small polaron transformation                                                                       
    max_occ_qu_vib_modes = [2]                                                            
    dim_el = 2**Nel # Number of fermionic Fock states
    dim_ph = np.prod(np.array(max_occ_qu_vib_modes,dtype=int)+1) # Number of bosonic Fock states
    dim_rho = dim_el*dim_ph # Number of states in the total fermionic-bosonic Fock space

    Single_El_Int_0 = np.array([[0.1]],dtype=float)                 # Energies of levels included in system, as well as hopping between levels
    Double_El_Int_0 = np.array([[0]],dtype=float)                   # Coulomb interactions between fermions - expressed as lower triangular matrix filled with two-particle interactions
    Vib_Freq_qu = np.array([0.02],dtype=float)                          # Vibrational frequency of phonon modes included in transport
    El_Ph_Int = np.array([[0.01]],dtype=float)                      # Electron-phonon coupling strength for the two types of couplings (bond length and site rigidity)   
    Single_El_Int = Single_El_Int_0 - (El_Ph_Int**2)/Vib_Freq_qu        # Fermionic energies after small polaron transformation

    if Nel > 1:
        Double_El_Int = Double_El_Int_0 - 2*(np.triu(np.matmul(El_Ph_Int.transpose(),El_Ph_Int)) + np.triu(np.matmul(El_Ph_Int.transpose(),El_Ph_Int)).transpose())/Vib_Freq_qu 
                                                                    # Coulomb fermionic interactions after small polaron transformation
    else:
        Double_El_Int = Double_El_Int_0
    
    Constraints = [Nel,N_qu_vib_modes,max_occ_qu_vib_modes]
    if N_qu_vib_modes == 0:  # Test whether the user wants to generate a fermionic or bosonic only Fock space, or a joint one.
        CreAnn1 = CreAnn.CreAnn(Constraints,'Fermi')
        d_ops,d,ddag,Fermionic_Fock_states = CreAnn1.return_operators()
    elif Nel == 0:
        CreAnn1 = CreAnn.CreAnn(Constraints,'Bose')
        a_ops,a,adag,Bosonic_Fock_states = CreAnn1.return_operators()
    elif N_qu_vib_modes != 0 and Nel != 0:
        CreAnn1 = CreAnn.CreAnn(Constraints,'Both') 
        d_ops,d,ddag,a_ops,a,adag,Both_Fock_states = CreAnn1.return_operators()

    Ham = np.zeros((dim_rho,dim_rho),dtype=float) # Initialize Hamiltonian (after small polaron transformation) to be filled
    Ham += Vib_Freq_qu[0]*np.matmul(adag[:,:,0],a[:,:,0]) # Fill Hamiltonian with bosonic energies (after small polaron transformation)
    for itrm1 in range(Nel): # Loop through fermionic levels 
        Ham += Single_El_Int[0,itrm1]*np.matmul(ddag[:,:,itrm1],d[:,:,itrm1])  # Add the energy of that fermionic level to the Hamiltonian
        for itrm2 in range(Nel): # Loop through all fermionic levels again to take into account double electron interactions
            Ham += Double_El_Int[itrm1,itrm2]*np.matmul(np.matmul(ddag[:,:,itrm1],d[:,:,itrm1]),np.matmul(ddag[:,:,itrm2],d[:,:,itrm2]))  # Include double electron interactions in Hamiltonian

    ham_file = open('Hamiltonian.txt',"w") # Open the Hamiltonian file with write access
    ham_file.write("-----------------------------------------------------------------------------------HAMILTONIAN-------------------------------------------------------------------------------\n")
    np.savetxt(ham_file,Ham,fmt='%3.2f') # Input the Hamiltonian
    ham_file.close()

    ph_file = open('Bosonic_Operators.txt',"w") # Do the same with all bosonic and fermionic annihilation and creation operators
    ph_file.write("-----------------------------------------------------------------------------------BOSONIC CREATION OPERATORS----------------------------------------------------------------------\n")
    np.savetxt(ph_file,adag[:,:,0],fmt='%4.2f')
    ph_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    np.savetxt(ph_file,adag[:,:,1],fmt='%4.2f')
    ph_file.write("---------------------------------------------------------------------------------BOSONIC ANNIHILATION OPERATORS--------------------------------------------------------------------\n")
    np.savetxt(ph_file,a[:,:,0],fmt='%4.2f')
    ph_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    np.savetxt(ph_file,a[:,:,1],fmt='%4.2f')
    ph_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    ph_file.close()

    el_file = open('Fermionic_Operators.txt',"w")
    el_file.write("-----------------------------------------------------------------------------------FERMIONIC CREATION OPERATORS----------------------------------------------------------------------\n")
    np.savetxt(el_file,ddag[:,:,0],fmt='%-2.1i')
    el_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    np.savetxt(el_file,ddag[:,:,1],fmt='%-2.1i')
    el_file.write("---------------------------------------------------------------------------------FERMIONIC ANNIHILATION OPERATORS--------------------------------------------------------------------\n")
    np.savetxt(el_file,d[:,:,0],fmt='%-2.1i')
    el_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    np.savetxt(el_file,d[:,:,1],fmt='%-2.1i')
    el_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    el_file.close()

######################################### CODE BELOW USES LANG-FIRSOV TRANSFORMATION ############################################################


    # Ham1 = np.zeros((dim_rho,dim_rho),dtype=float) 
    # Ham1 += Single_El_Int[0,0]*np.matmul(ddag[:,:,0],d[:,:,0]) + Vib_Freq_qu[0]*np.matmul(adag[:,:,0],a[:,:,0])

    # Ham = np.zeros((dim_rho,dim_rho),dtype=float) # Initialize Hamiltonian (after small polaron transformation) to be filled
    # Ham += Vib_Freq_qu[0]*np.matmul(adag[:,:,0],a[:,:,0]) # Fill Hamiltonian with bosonic energies (after small polaron transformation)
    # for itrm1 in range(Nel): # Loop through fermionic levels 
    #     Ham += Single_El_Int[0,itrm1]*np.matmul(ddag[:,:,itrm1],d[:,:,itrm1])  # Add the energy of that fermionic level to the Hamiltonian
    #     for itrm2 in range(Nel): # Loop through all fermionic levels again to take into account double electron interactions
    #         Ham += Double_El_Int[itrm1,itrm2]*np.matmul(np.matmul(ddag[:,:,itrm1],d[:,:,itrm1]),np.matmul(ddag[:,:,itrm2],d[:,:,itrm2]))  # Include double electron interactions in Hamiltonian

    # # Generate fermionic creation and annihilation operators dressed with Franck-Condon matrix from Franck_Condon class

    # FC_Operators = Franck_Condon.Franck_Condon(Constraints,El_Ph_Int,Vib_Freq_qu) # Run Franck_Condon program to generate Franck-Condon matrix and the fermionic creation and annihilation operators
    # #                                                                           # dressed with the Franck-Condon matrix. 
    # FC_Matrix,FC_Matrix_Fock_Space = FC_Operators.return_FC_Operators()

    # d_FC = np.zeros((dim_rho,dim_rho,Nel))
    # ddag_FC = np.zeros((dim_rho,dim_rho,Nel))
    # d_ops_FC = np.zeros((dim_rho,dim_rho,Nel,2))
    # for itrm in range(Nel):
    #     d_FC[:,:,itrm] = np.matmul(d[:,:,itrm],FC_Matrix_Fock_Space[:,:,itrm])
    #     ddag_FC[:,:,itrm] = np.matmul(ddag[:,:,itrm],np.transpose(FC_Matrix_Fock_Space[:,:,itrm]))
    #     d_ops_FC[:,:,itrm,0] = d_FC[:,:,itrm]
    #     d_ops_FC[:,:,itrm,1] = ddag_FC[:,:,itrm]
    
    # return d_ops_FC,d_FC,ddag_FC,FC_Matrix,a_ops,a,adag,Both_Fock_states,Ham 
    # return d_ops,d,ddag,a_ops,a,adag,Both_Fock_states,Ham 


###################### JUNK - THIS IS DONE IN THE CreAnn FUNCTION NOW ##############################

# def fock_states(Single_El_Int,Vib_Freq_qu,Nel,N_qu_vib_modes,max_occ_qu_vib_modes,dim_ph,dim_el):

#     # Generate associated Fock states using basis of electron and phonon occupancy
#     """
#     This ordering follows that of the fermionic and bosonic annihilation and creation operators
#     """

#     El_Occ = np.arange(0,2); El_Occ.astype(int); El_Occ.shape = (2,1)
#     El_States_1 = np.tile(np.tile(El_Occ,(Nel,1)),(dim_ph,1))
#     # El_States_2 = np.tile(np.matrix.repeat(El_Occ,Nel,axis=0),(dim_ph,1))
#     Ph_States_1 = np.arange(0,max_occ_qu_vib_modes[0]+1); Ph_States_1.astype(int); Ph_States_1.shape = (max_occ_qu_vib_modes[0]+1,1)
#     Ph_States_1 = np.matrix.repeat(Ph_States_1,dim_el,axis=0)
#     # Ph_States_2 = np.arange(0,max_occ_qu_vib_modes[1]+1); Ph_States_2.astype(int); Ph_States_2.shape = (max_occ_qu_vib_modes[1]+1,1)
#     # Ph_States_2 = np.tile(np.matrix.repeat(Ph_States_2,dim_el,axis=0),(max_occ_qu_vib_modes[0]+1,1))

#     # ElPh_Fock_States = np.concatenate((El_States_1,El_States_2,Ph_States_1,Ph_States_2),axis=1) 
#     ElPh_Fock_States = np.concatenate((El_States_1,Ph_States_1),axis=1) 

#     Energy_El_1 = Single_El_Int[0,0]*El_States_1
#     Energy_El_2 = Single_El_Int[1,1]*El_States_2
#     Energy_Ph_1 = Vib_Freq_qu[0]*Ph_States_1
#     Energy_Ph_2 = Vib_Freq_qu[1]*Ph_States_2
#     diag_Energies = np.sum(np.concatenate((Energy_El_1,Energy_El_2,Energy_Ph_1,Energy_Ph_2),axis=1),axis=1)
    

