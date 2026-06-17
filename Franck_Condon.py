import numpy as np
import itertools
import scipy.special
from scipy.special import comb
import CreAnn
import matplotlib.pyplot as plt

class Franck_Condon():
    
    def __init__(self,Constraints,El_Ph_Int,Ph_Freq):
        
        self.Constraints = Constraints
        self.Nel = self.Constraints[0]
        self.max_ph_occ = self.Constraints[2][0]
        self.El_Ph_Int = El_Ph_Int
        self.Ph_Freq = Ph_Freq
        
        self.Laguerre_Polynomials()
        self.Franck_Condon_Matrix()
        
    def Laguerre_Polynomials(self):
        x = ((self.El_Ph_Int)**2)/((self.Ph_Freq)**2)
        self.Laguerre = np.zeros((self.max_ph_occ+1,self.max_ph_occ+1,self.Nel))
        for itrm in range(self.Nel):
            self.Laguerre[0,:,itrm] = 1
            self.Laguerre[1,:,itrm] = -x[itrm] + np.arange(self.max_ph_occ+1) + 1
            for n in range(2,self.max_ph_occ+1):
                self.Laguerre[n,:,itrm] = ((2*(n-1) + np.arange(self.max_ph_occ+1) + 1 - x[itrm])*self.Laguerre[n-1,:,itrm] - (n-1+np.arange(self.max_ph_occ+1))*self.Laguerre[n-2,:,itrm])/n

    def Franck_Condon_Matrix(self):
        dim_ph = self.max_ph_occ+1
        dim_el = 2**(self.Nel)
        dim_rho = dim_el*dim_ph
        self.FC_Matrix = np.zeros((dim_ph,dim_ph,self.Nel))
        self.FC_Matrix_Fock_Space = np.zeros((dim_rho,dim_rho,self.Nel))
        for itrm in range(self.Nel):
            for nurow in range(self.max_ph_occ+1):
                for nucol in range(self.max_ph_occ+1):
                    self.FC_Matrix[nurow,nucol,itrm] = np.exp(-0.5*((self.El_Ph_Int[itrm]/self.Ph_Freq)**2))*np.sqrt(np.math.factorial(np.amin([nurow,nucol]))/np.math.factorial(np.amax([nurow,nucol])))\
                                            *((np.sign(nucol-nurow)*(self.El_Ph_Int[itrm]/self.Ph_Freq))**(np.abs(nucol-nurow)))*self.Laguerre[np.amin([nurow,nucol]),np.abs(nucol-nurow),itrm]
            self.FC_Matrix_Fock_Space[:,:,itrm] = np.kron(self.FC_Matrix[:,:,itrm].transpose(),np.eye(dim_el))

    def return_FC_Operators(self):
        # return self.D_ops_FC,self.d_FC,self.ddag_FC,self.FC_Matrix
        return self.FC_Matrix,self.FC_Matrix_Fock_Space 

    def Dressed_FD_Functions(self,e,w,T):
        q = [np.arange(self.max_ph_occ+1)]
        FC_Sq = self.FC_Matrix[:,:,0]**2
        self.FD_01 = FC_Sq*(1 - 1/(1+np.exp((e - w*(np.transpose(q)-q))/T)))
        self.FD_10 = FC_Sq*(1/(1+np.exp((e + w*(np.transpose(q)-q))/T)))
        return self.FD_01,self.FD_10

if __name__=='__main__':
    
    import Franck_Condon

    Nel = 1
    Nph = 1
    max_ph_occ = 5    
    Constraints = np.array([Nel,Nph,[max_ph_occ]])
    El_Ph_Int = np.array([[1,3]])
    Ph_Freq = 1

    FC_Operators = Franck_Condon.Franck_Condon(Constraints,El_Ph_Int,Ph_Freq)
    # D_ops_FC,d_FC,ddag_FC,FC_Matrix = FC_Operators.return_FC_Operators()
    FC_Matrix,FC_Matrix_Fock_Space = FC_Operators.return_FC_Operators()

    FC_Operators_file = open('Franck_Condon_Operators.txt',"w")

    # FC_Operators_file.write("-----------------------------------------------------------------------------------FERMIONIC CREATION OPERATORS----------------------------------------------------------------------\n")
    # for itrm in range(Constraints[0]):
    #     np.savetxt(FC_Operators_file,ddag_FC[:,:,itrm],fmt='%-5.5f')
    #     FC_Operators_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

    # FC_Operators_file.write("-----------------------------------------------------------------------------------FERMIONIC ANNIHILATION OPERATORS----------------------------------------------------------------------\n")
    # for itrm in range(Constraints[0]):
    #     np.savetxt(FC_Operators_file,d_FC[:,:,itrm],fmt='%-5.5f')
    #     FC_Operators_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

    FC_Operators_file.write("-----------------------------------------------------------------------------------FRANCK-CONDON MATRICES----------------------------------------------------------------------\n")
    for itrm in range(Constraints[0]):
        np.savetxt(FC_Operators_file,FC_Matrix[:,:,itrm],fmt='%-5.5f')
        FC_Operators_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

    FC_Operators_file.close()

    FC_Operators_Fock_Space_file = open('Franck_Condon_Operators_Fock_Space.txt',"w")

    FC_Operators_Fock_Space_file.write("-----------------------------------------------------------------------------------FRANCK-CONDON MATRICES----------------------------------------------------------------------\n")
    for itrm in range(Constraints[0]):
        np.savetxt(FC_Operators_Fock_Space_file,FC_Matrix_Fock_Space[:,:,itrm],fmt='%-5.5f')
        FC_Operators_Fock_Space_file.write("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

    FC_Operators_Fock_Space_file.close()

    for itrm in range(Constraints[0]):
        plt.figure()
        FC_Matrix_Sq = (FC_Matrix[:,:,itrm])**2
        c = plt.imshow(FC_Matrix_Sq,cmap='inferno',extent=[0,max_ph_occ,0,max_ph_occ],origin='lower')
        plt.colorbar(c)

    plt.show()