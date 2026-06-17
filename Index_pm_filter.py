# ------------------------------------------------------------------------------------
#
#                   GENERATING INDICES OF THE HQME
#
# ------------------------------------------------------------------------------------
#
# This class contains the functions necessary to generate the unique indices of each tier 
# of the HQME. It also generates arrays that allow the HQME to quickly find ADOs
# connected between the tiers. Finally, it has functionality to include the importance criterion
# defined in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.235426 . 
# 
# Samuel Rudge
#
# USAGE: Indices = Index_pm_filter.Hierarchy_index(Nmax,Nel,Lmax,Nleads,Nsign,Nmodes)
#
# INPUTS:
#       Nmax                                    Maximum tier of the hierarchy
#
#       Nel                                     Number of electronic levels
#
#       Lmax                                    Number of Pade poles included in the Fermi-Dirac approximation
#
#       Nleads                                  Number of leads system is attached to; usually 2 for fermionic baths
#
#       Nsign                                   Refers to \sigma = +/-, which is shorthand for creation/annihilation operators.
#                                                Always 2.
#
#       Nmodes                                  Total number of j = {K,\sigma,l,m indices} indices. (Called modes for convention, not related
#                                               to vibrational modes. Calculated as Nel*(Lmax+1)*Nleads*Nsign)
#
# OUTPUT:
#       Indices                                 Object containing functions to return the index information and to apply the importance criterion,
#                                               if desired.

import numpy as np                                                              # Import necessary functions
import itertools
import sys
import copy
# import numba as nb
# from numba import int64, int32    # import the types
# from numba.experimental import jitclass

# spec = [('Nmax',int32),('Lmax',int32),('Nel',int32),('Nleads',int32),('Nsign',int32),('Npoles',int32),('Nmodes',int32),('Un_Ind',int32[:,:]),
#         ('Hier_ind',int32[:,:,]),('KsigLm',int32[:,:,]),('Old',int32[:,:])]

# @jitclass(spec)
class Hierarchy_index():                                                                    

    def __init__(self,Nmax,Nel,Lmax,Nleads,Nsign,Nmodes,wbl_YN):                       # Run necessary functions when 'Indices' object is created.
       
        self.Nmax = int(Nmax)                                                   # Define parameters from inputs
        self.Lmax = int(Lmax)
        self.Nel = int(Nel)
        self.Nleads = int(Nleads)
        self.Nsign = int(Nsign)
        self.wbl_YN = wbl_YN
        if self.wbl_YN == 1:
            self.Npoles = int(self.Lmax)                                        # Total number of Pade poles is Lmax + 1
        elif self.wbl_YN == 0:
            self.Npoles = int(self.Lmax + 1)                                        # Total number of Pade poles is Lmax + 1
        else:
            raise ValueError("spectral density can only be wide-band limit (wbl_YN == 'TRUE') or not wide-band limit (wbl_YN == 'FALSE')")
        self.Nmodes = int(Nmodes)
        
        self.Gen_Un_Ind()                                                       # Run Gen_Un_Ind(),Gen_Ind_Plus(),Gen_Ind_Minus(),Gen_KsigLm() automatically when object is created
        self.Gen_Ind_Plus()
        self.Gen_Ind_Minus()
        self.Gen_KsigLm()

    # --------------------------------------------------------------------- 
    #
    #                     GENERATE ALL INDEX MODES
    #
    # ---------------------------------------------------------------------
    #
    # This function generates an array of size [Nmodes,1] containing all index modes: {[0:Nleads,0:Nsig,0:Lmax,0:Nel]}.
    # For example, for the parameters Nleads=2 (left=0,right=1); Nsign=2(+=0,-=1); Lmax=1; Nel=1, the index modes are 
    # modes: 0:[0,0,0,0]    1:[0,1,0,0]    2:[0,0,1,0]    3:[0,1,1,0]      
    #        4:[1,0,0,0]    5:[1,1,0,0]    6:[1,0,1,0]    7:[1,1,1,0]      
    # Nmodes=8
    #
    # USAGE: KsigLm = Indices.Gen_KsigLm()
    #
    # INPUTS:
    #        No inputs required that have not already been provided when creating 'Indices' object
    #
    # OUTPUT: 
    #        KsigLm                             Array containing index values of K,sigma,l,m for each mode (see description for example)
    
    def Gen_KsigLm(self):                                                       
    
        Mu = np.arange(0,self.Nleads); Mu.astype(int); Mu.shape = (self.Nleads,1); Mu = np.matrix.repeat(Mu,self.Nsign*self.Npoles*self.Nel,axis=0) 
        Sig = np.arange(0,self.Nsign); Sig.astype(int); Sig.shape = (self.Nsign,1); Sig = np.tile(Sig,(self.Nleads*self.Npoles*self.Nel,1))
        Poles = np.arange(0,self.Npoles); Poles.astype(int); Poles.shape = (self.Npoles,1); Poles = np.tile(np.matrix.repeat(Poles,self.Nsign,axis=0),(self.Nleads*self.Nel,1))
        Electronic_States = np.arange(0,self.Nel); Electronic_States.astype(int); Electronic_States.shape = (self.Nel,1)
        Electronic_States = np.tile(np.matrix.repeat(Electronic_States,self.Nsign*self.Npoles,axis=0),(self.Nleads,1))
        self.KsigLm = np.concatenate((Mu,Sig,Poles,Electronic_States),axis=1); self.KsigLm.astype(int) 
                                                                                # Concatenate four [Nmodes,1] arrays into one array. 

        if self.wbl_YN == 1:
            Mu = np.arange(0,self.Nleads); Mu.astype(int); Mu.shape = (self.Nleads,1); Mu = np.matrix.repeat(Mu,self.Nsign*self.Nel,axis=0) 
            Sig = np.arange(0,self.Nsign); Sig.astype(int); Sig.shape = (self.Nsign,1); Sig = np.tile(Sig,(self.Nleads*self.Nel,1))
            Electronic_States = np.arange(0,self.Nel); Electronic_States.astype(int); Electronic_States.shape = (self.Nel,1)
            Electronic_States = np.tile(np.matrix.repeat(Electronic_States,self.Nsign,axis=0),(self.Nleads,1))
            self.Ksig0m = np.concatenate((Mu,Sig,Electronic_States),axis=1); self.KsigLm.astype(int) 

    # ---------------------------------------------------------------------
    #
    #                  GENERATE UNIQUE INDICES FOR EACH TIER
    #
    # ---------------------------------------------------------------------
    #
    # This function generates indices corresponding to the unique and non-zero ADOs, after taking into account 
    # the hermiticity relation and the multiple application of the same Grassmann number in one ADO. See the notes 
    # for a detailed explanation of why each ADO is retained or removed.
    #
    # USAGE: Un_Ind,Hier_ind,tier_ind,len_un_index,len_index_plus = Indices.Gen_Un_Ind()
    #
    # INPUTS: 
    #        No inputs required that have not already been provided when creating 'Indices' object
    #
    # OUTPUTS:
    #        len_un_index                        Number of unique ADOs 
    #
    #        Un_Ind                              Array containing the modes used to build each ADO. It has size [len_un_index,Nmax].
    #                                            Each row corresponds to one of the unique ADOs used in the hierarchy, arranged in 
    #                                            order of increasing tier (e.g. the zeroth tier is listed 1st, then the 1st tier,
    #                                            then the 2nd, and so on.) Each column corresponds to the modes used to create that 
    #                                            ADO, in reverse order. If an ADO is of a lower tier than Nmax, its remaining columns 
    #                                            are filled with -1. For example, if Nmax = 4, then the 1st row, corresponding to the 
    #                                            zeroth tier, is [-1,-1,-1,-1]. If we further assume that all parameters are the same as 
    #                                            in the example in Gen_KsigLm, then the next four rows are the 1st tier ADOs:
    #                                            [0,-1,-1,-1],[2,-1,-1,-1],[4,-1,-1,-1],[6,-1,-1,-1]. 
    #                                            The 2nd tier ADOs are 
    #                                            [0,1,-1,-1],[0,2,-1,-1],[0,3,-1,-1],...,[0,7,-1,-1]
    #                                            [2,3,-1,-1],[2,4,-1,-1],...,[2,7,-1,-1]
    #                                            ...
    #                                            [6,7,-1,-1]
    #                                            The third tier ADOs are 
    #                                            [0,1,2,-1],[0,1,3,-1],[0,1,4,-1],...,[0,1,7,-1]
    #                                            [0,2,3,-1],[0,2,4,-1],...,[0,2,7,-1]
    #                                            ...
    #                                            [2,3,4,-1],[2,3,5,-1],...,[2,3,7,-1]
    #                                            [2,4,5,-1],[2,4,6,-1],[2,4,7,-1].
    #                                            ...
    #                                            [5,6,7,-1]
    #                                            And the fourth tier ADOs are
    #                                            [0,1,2,3],[0,1,2,4],[0,1,2,5],[0,1,2,6],[0,1,2,7]
    #                                            ...
    #                                            [4,5,6,7].
    #
    #       Hier_ind                             Array with size [Nmax+1,2] containing the start and end index of each tier. For example,
    #                                            the zeroth tier starts and ends at index 0, so the 1st row is always [0,0]. The 2nd
    #                                            row corresponds to the start and end index of the 1st tier ADOs, which must be [1,Nmodes/2].
    #                                            This continues up to the Nmax-th tier (row).
    #
    #       tier_ind                             Array with size [1,len_un_ind] recording the tier of each ADO listed in Un_Ind
    #
    #       len_index_plus                       The total number of ADOs that connect to the tier above (Basically all ADOs except the maximum tier). Used in 
    #                                            the Gen_Ind_Plus function.

    def Gen_Un_Ind(self):
                
        self.Un_Ind = np.full((1,self.Nmax),-1,dtype=int)                       # Initialize Un_Ind matrix with more rows than necessray (at this point we do not know
                                                                                    # how many rows (ADOs) it needs to have). It is currently filled with -1
        Un_Ind_tier_1 = np.full((int(self.Nmodes/2),self.Nmax),-1,dtype=int)
        Un_Ind_tier_1[:,0] = np.arange(0,self.Nmodes-1,2)             # Fill the 1st tier ADO modes in
        self.Un_Ind = np.append(self.Un_Ind,Un_Ind_tier_1,axis=0)

        self.Hier_ind = np.empty((self.Nmax+1,2),dtype=int)                         # Initialize Hier_Ind matrix with zeros
        self.Hier_ind[0,0:2:1] = [int(0),int(0)]                                    # Fill the 1st row in with the indices of the zeroth tier (always [0,0])
        self.Hier_ind[1,0:2:1] = [int(1),int((self.Nmodes//2))]                     # Fill the 2nd row in with the indices of the 1st tier
        indcount = int(1 + (self.Nmodes//2))                                        # Create variable indcount to keep track of number of unique ADOs

        for itrn in range(2,self.Nmax+1):                                           # Loop through tiers 2 and up (0 and 1 already filled in)
            for itrnj in range(self.Hier_ind[itrn-1,0],self.Hier_ind[itrn-1,1]+1):  # Loop through indices of the previous tier to build the unique ADOs of the current tier 
                Old = np.repeat(np.array([self.Un_Ind[itrnj,0:(itrn-1):1]]),int(self.Nmodes-1-self.Un_Ind[itrnj,itrn-2]),0) 
                                                                                    # Build an array containing repeats of the modes in the ADO of the previous tier with 
                                                                                    # index itrnj. The number of repeats is just the number of modes greater than the last 
                                                                                    # mode in the list up to and including the maximum mode number. E.g. if itrn = 2 and 
                                                                                    # itrnj = 3 (assuming the same parameters as previously), then Old = [[4],[4],[4]], corresponding
                                                                                    # to the 2nd tier ADOs [4,5],[4,6],[4,7].
                New = np.arange(self.Un_Ind[itrnj,itrn-2]+1,self.Nmodes); New.shape = ((int(self.Nmodes-1-self.Un_Ind[itrnj,itrn-2]),1))
                                                                                    # Create a column vector containing the modes added to the previous tier to create the new tier.
                                                                                    # For the previous example, New = [[4],[5],[6]]
                array_to_append = np.full((int(self.Nmodes-1-self.Un_Ind[itrnj,itrn-2]),self.Nmax),-1,dtype=int)
                Old_New_Combined = np.concatenate((Old,New),axis=1)
                array_to_append[:,0:itrn] = Old_New_Combined
                self.Un_Ind = np.append(self.Un_Ind,array_to_append,axis=0)
                                                                                    # Concatenate these two and add them to the Un_Ind array (e.g. [[4,5],[4,6],[4,7]])
                indcount = int(indcount + (self.Nmodes - 1 - self.Un_Ind[itrnj,(itrn-2)])) 
                                                                                    # Calculate the number of new ADOs generated by this procedure and update the total number
            self.Hier_ind[itrn,0] = self.Hier_ind[itrn-1,1]+1; self.Hier_ind[itrn,1] = indcount-1
                                                                                    # Include the indices of the new ADOs of this tier in Hier_Ind.

        self.len_un_ind = indcount                                                  # Define the total number of unique ADOs (= unique indices)
        self.len_index_plus = self.Hier_ind[self.Nmax][0]                           # The total number of ADOs that connect to a tier above (exclude the Nmax-th tier)
        self.tier_index = np.empty(self.len_un_ind,dtype=int)                       # Initialize array of tier indices
        for itrn in range(0,self.Nmax+1):                                           # Loop through the tiers
            tier_index_lower = self.Hier_ind[itrn,0]                                # Find the lower and upper indices of that tier
            tier_index_upper = self.Hier_ind[itrn,1]+1
            self.tier_index[tier_index_lower:tier_index_upper] = itrn               # Fill the tier_index array

    # ---------------------------------------------------------------------
    #
    #          CREATE ARRAY CONNECTING nth TIER TO (n+1)th TIER
    #
    # ---------------------------------------------------------------------
    #
    # This function generates an array that details the connection between ADOs of tier n 
    # to ADOs of tier n+1, allowing for rapid indexing in the HQME. Consider the general problem:
    # You have an ADO of tier n, defined by modes {j_{n},...,j_{1}}. You know that the HQME connects it 
    # to ADOs of tier n+1, defined by modes {j_{n+1},j_{n},...,j_{1}}. But to include this connection 
    # numerically, you must know (a) the index of this ADO and (b) the index of the ADO it connects to.
    # You could run a search on Un_Ind to find the corresponding row, but to do this for all ADOs and all
    # ADO connections would be expensive and slow. The idea of the array that this function creates is that 
    # it takes the modes defining the nth tier ADO and the mode you want to add to connect it to the (n+1)th 
    # tier ADO, and it quickly can find the corresponding indices.
    #
    # USAGE: Index_plus = Indices.Gen_Ind_Plus()
    #
    # INPUTS: 
    #        No inputs required that have not already been provided when creating 'Indices' object
    #
    # OUTPUTS:
    #        Index_plus                         Array of size [len_ind_plus,Nmodes]. Each row corresponds to an ADO and each column represents
    #                                           the Nmodes possible ADOs that can be built from this ADO; i.e. all ADOs of a higher tier
    #                                           that connect to it in the HQME. Each element is actually an array of size [,3]. The 1st element
    #                                           contains the index of the UNIQUE ADO it connects to. The 2nd element contains the number of permutations
    #                                           required to connect the two, and the third element records whether you need to take the hermitian conjugate 
    #                                           to connect them (0 = no hermitian conjugate and 1 = hermitian conjugate). These latter two elements are 
    #                                           necessary to keep track of the (-1) prefactors that must be taken into account when removing the linear 
    #                                           dependencies in the HQME. As an example, consider the same example as before. The zeroth tier connects to 
    #                                           the 1st tier ADOs, 
    #                                           so the elements of the 1st row are 
    #                                           [1,0,0] = You can find the 1st tier ADO defined by j_{0} by going to ADO index 1 (Row 1 in Un_Ind)
    #                                           [1,0,1] = You can find the 1st tier ADO defined by j_{1} by going to ADO index 1 and using the hermiticity relationship
    #                                           [2,0,0] = You can find the 1st tier ADO defined by j_{2} by going to ADO index 2 (Row 2 in Un_Ind)
    #                                           [2,0,1] = You can find the 1st tier ADO defined by j_{3} by going to ADO index 2 and using the hermiticity relationship
    #                                           ...
    #                                           [4,0,0] = You can find the 1st tier ADO defined by j_{6} by going to ADO index 4 (Row 4 in Un_Ind)
    #                                           [4,0,1] = You can find the 1st tier ADO defined by j_{7} by going to ADO index 4 and using the hermiticity relationship
    #                                           The 2nd row elements (corresponding to the 1st tier ADO with mode j_{0}) are 
    #                                           [-1,-1,-1] = The 2nd tier ADO with modes {j_{0},j_{0}} can be excluded from the hierarchy
    #                                           [5,0,1] = You can find the 2nd tier ADO defined by {j_{1},j_{0} by going to ADO index 5 (Row 5 in Un_Ind)
    #                                           ...
    #                                           [11,0,1] = You can find the 2nd tier ADO defined by {j_{1},j_{0} by going to ADO index 11 (Row 11 in Un_Ind)
    #                                           The third row elements (corresponding to the 1st tier ADO with mode j_{2}) are
    #                                           [6,1,0] = You can find the 2nd tier ADO defined by {j_{0},j_{2}} by going to ADO index 6 and permuting once
    #                                           [7,1,1] = You can find the 2nd tier ADO defined by {j_{1},j_{2}} by going to ADO index 7, permuting once, and applying the 
    #                                            hermiticity relationship
    #                                            ...
    #                                            The 20th row elements (corresponding to the 2nd tier ADO with modes {j_{7},j_{6}}) are 
    #                                           [41,2,0] = You can find the 3rd tier ADO with modes {j_{0},j_{7},j_{6}} by going to ADO index 41 (Row 41 in Un_Ind) and 
    #                                           permuting twice
    #                                           ...
    #                                           [54,3,1] = You can find the 3rd tier ADO with modes {j_{5},j_{7},j_{6}} by going to ADO index 54 (Row 54 in Un_Ind) and 
    #                                            permuting three times and applying the hermiticity relation
    #                                           [-1,-1,-1] = You can exclude the 3rd tier ADO with modes {j_{6},j_{7},j_{6}} from the hierarchy
    #                                           [-1,-1,-1] = You can exclude the 3rd tier ADO with modes {j_{7},j_{7},j_{6}} from the hierarchy

    def Gen_Ind_Plus(self):
        
        self.Index_plus = [[[]]]*(self.Hier_ind[self.Nmax-1,1]+1)                   # Initialize Index_plus as empty
        Index_plus_0 = [[]*3]*self.Nmodes                                           # Initialize the elements of the first row (for zeroth tier) as empty
        plus_ind = 0                                                                # Initialize varible that keeps track of which unique (???)

        for itrj0 in range(self.Nmodes):                                            # Loop through all modes to generate elements of first row
            plus_list = [0,0,0]                                                     # Put the first element in Index_plus_0 in
            if itrj0 % 2 != 0:                                                      # Check if this mode is odd
                plus_list[0] = plus_ind                                             # If mode is odd, make ADO index equal to current value of plus_ind
                plus_list[2] = 1                                                    # If mode is odd, make hermiticity element equal to 1
                Index_plus_0[itrj0] = plus_list                                     # Fill the itrj0-th column with this array: [plus_ind,0,1]
            elif itrj0 % 2 == 0:                                                    # Check if this mode is even
                plus_ind += 1                                                       # Add 1 to the plus_ind
                plus_list[0] = plus_ind                                             # Make ADO index equal to current value of plus_ind
                Index_plus_0[itrj0] = plus_list                                     # Fill the itrj0-th colum with this array: [plus_ind,0,1]
             
        self.Index_plus[0] = Index_plus_0                                           # Fill the first column of Index_plus with Index_plus_0

        for itrn in range(1,self.Nmax):                                             # Loop through all tiers except Nmax
            for itrjn in range(self.Hier_ind[itrn,0],self.Hier_ind[itrn,1]+1):      # Loop through all ADOs of this tier (itrjn is the index of the corresponding ADO)
                j_old = self.Un_Ind[itrjn,0:itrn:1]                                 # Create vector of modes corresponding to the ADO of this index
                Index_plus_n = [[]]*self.Nmodes                                     # Initialize the elements of the itrjn-th row of Index_plus as empty
                for itrjnp1 in range(self.Nmodes):                                  # Loop through all modes; they will be appended to this ADO to make the ADO of the (n+1)th tier
                    j_new = np.append(j_old,itrjnp1)                                # Create vector of modes corresponding to ADO of tier n+1
                    if itrjnp1 in j_old:                                            # Check if this ADO contains the same Grassmann number twice
                        Index_plus_n[itrjnp1] = [-1,-1,-1]                          # If it does, exclude it from hierarchy by putting -1s in Index_plus
                    elif itrjnp1 < j_old[itrn-1]:                                   # Check if the new mode added to the ADO is less than the last mode in j_old
                        j_same = np.sort(j_new)                                     # If it is, it can be rewritten using permutations and the hermiticity relation
                        nperm = len(np.extract(j_old > itrjnp1,j_old))              # Calculate number of permutations required to put j_new in ascending order
                        if j_same[0] % 2 == 0:                                      # Check if new mode is even 
                            ind_jsame = 0                                           # If mode is even, there is no need to apply the hermiticity relation, so initialize index of permuted ADO
                            for itrjsame in range(itrn+1):                          # Loop through modes of permuted ADO in order to find index
                                ind_jsame = self.Index_plus[ind_jsame][int(j_same[itrjsame])][0]
                                                                                    # Use fast indexing with already completed part of Index_plus (by defn it must be already completed)
                                                                                    # e.g. if the permuted ADO has modes {j_{3},j_{2},j_{0}} this part of the code goes first to the zeroth
                                                                                    # row of Index_plus (corresponding to the zeroth tier) and searches for the 0th column, which is the element
                                                                                    # containing the connection to the first tier ADO with mode {j_{0}}. The zeroth element of this element
                                                                                    # is the index of this first tier ADO, which becomes the new ind_jsame. Then, the code goes to the row of 
                                                                                    # Index_plus corresponding to this new first tier ADO, and the column corresponding to j_{2}; thus finding the 
                                                                                    # index of {j_{2},j_{0}}. Then it repeats this to find the index of {j_{3},j_{2},j_{0}}, systematically 
                                                                                    # updating the index until it has the correct one for jsame. This is much faster than running a search on 
                                                                                    # Un_Ind every time.
                            plus_list = [ind_jsame,nperm,0]                         # The old ADO connects to the new ADO via nperm permutations of the ADO with index ind_jsame. There are no 
                                                                                    # hermiticity relations necessary, so that element remains 0.
                            Index_plus_n[itrjnp1] = plus_list                       # This new information is inserted into the list of connections for this ADO.
                        elif j_same[0] % 2 != 0:                                    # If the newly added mode is odd, then we need to do the same steps but add the hermiticity relation
                            ind_evodd = np.arange(itrn+1,dtype=int)                 # Generate array of integers 0:n the same length as j_same (they will serve as indices for j_same)
                            ind_ev = ind_evodd[j_same%2==0]                         # Create array of indices of j_same corresponding to even modes
                            ind_odd = ind_evodd[j_same%2!=0]                        # Create array of indices of j_same corresponding to odd modes
                            j_same[ind_ev] += 1; j_same[ind_odd] -= 1               # Take Hermitian conjugate of j_same (odd indices go down by 1 and even indices go up by 1)
                            for itrpair in range(0,itrn):                           # Conjugate j_same may have modes out of order, so need to resort. Loop through each pair of elements
                                if j_same[itrpair+1] < j_same[itrpair]:             # If the pair of elements is out of order, switch them and add 1 to nperm
                                    j_same[[itrpair,itrpair+1]] = j_same[[itrpair+1,itrpair]]
                                    nperm += 1
                            ind_jsame = 0                                           # Finally have the unique ADO that j_old connects to. Need to do the same process to calculate its index.
                            for itrjsame in range(itrn+1): 
                                ind_jsame = self.Index_plus[ind_jsame][int(j_same[itrjsame])][0]
                            plus_list = [ind_jsame,nperm,1]                         # Construct element containing index of connecting ADO, number of permutations required, and 1 to indicate
                                                                                    # we had to perform a hermiticity relation
                            Index_plus_n[itrjnp1] = plus_list                       # This new information is inserted into the list of connections for this ADO.
                    else:
                        plus_ind += 1                                               # If the new ADO formed by j_new from j_old is unique (i.e. no permutations or hermiticity relation required),
                                                                                    # then by construction it will be the next unique ADO in Un_Ind, so add 1 to plus_ind.
                        plus_list = [plus_ind,0,0]                                  # Construct element containing index of connecting ADO with zero permutations and no hermiticity relation
                        Index_plus_n[itrjnp1] = plus_list                           # Insert new information into list of connection for this ADO
                self.Index_plus[itrjn] = Index_plus_n                               # Insert row of connection information for this ADO into the appropriate row of Index_plus

    # ---------------------------------------------------------------------
    #
    #          CREATE ARRAY CONNECTING nth TIER TO (n-1)th TIER
    #
    # ---------------------------------------------------------------------
    #
    # This function creates a similar array to Gen_Ind_Plus, except that now it details connectiosn between
    # ADOs of tier n and ADOs of tier n-1.
    #
    # USAGE: Index_minus = Indices.Gen_Ind_Minus()
    #
    # INPUTS: 
    #        No inputs required that have not already been provided when creating 'Indices' object
    #
    # OUTPUTS:
    #        Index_minus                        Array of size [len_un_ind,Nmax]. Each row corresponds to an ADO and each column represents
    #                                           the possible ADOs that can be built by removing a mode from this ADO; i.e. all ADOs of a lower tier
    #                                           that connect to it in the HQME. Evidently an ADO of the nth tier only fills up to column n < Nmax.
    #                                           Each element is actually an array of size [,4]. The 1st element contains the mode that is removed to make the connecting ADO.
    #                                           The 2nd element contains the index of the UNIQUE ADO it connects to. The 3rd element contains the number of permutations
    #                                           required to connect the two, and the 4th element records whether you need to take the hermitian conjugate 
    #                                           to connect them (0 = no hermitian conjugate and 1 = hermitian conjugate). Consider the same example as before. 
    #                                           The zeroth tier connects to no ADOs in a lower tier (there are no lower tiers) so all elements in the first row are 
    #                                           [-1,-1,-1,-1].
    #                                           The 1st tier ADOs connect to the zeroth tier by removing the only mode they contain, so the elements of the next few rows are  
    #                                           [0,0,0,0],[-1,-1,-1,-1],[-1,-1,-1,-1] = You can find the 0th tier ADO by removing mode j_{0} from this ADO
    #                                           and going to ADO index 0 (Row 0 in Un_Ind). Since it only has 1 mode, the other two elements are filled with -1.
    #                                           [2,0,0,0],[-1,-1,-1,-1],[-1,-1,-1,-1] = You can find the 0th tier ADO by removing mode j_{2} from this ADO
    #                                           and going to ADO index 0 (Row 0 in Un_Ind). Since it only has 1 mode, the other two elements are filled with -1.
    #                                           ...
    #                                           (Row 5) [0,1,0,1],[1,1,0,0],[-1,-1,-1,-1] = The row corresponds to a 2nd tier ADO with ADO index 5 (if you check 
    #                                           in Un_Ind this is defined by modes {j_{1},j_{0}}). In the first element, [0,1,0,1], the mode j_{0} is removed,
    #                                           (there is a 0 in the first place) this generates a 1st tier ADO with remaining mode 
    #                                           j_{1}, which can be expressed as the 1st tier ADO with mode j_{0} and ADO index 1 (the second element). Since we need to do a hermiticity 
    #                                           relation, a 1 goes in the last place. The second element, [1,1,0,0], corresponds to the same 2nd tier ADO with index 5, but now with 
    #                                           the j_{1} mode removed (there is a 1 in the first place) then you connect to the 1st tier ADO with ADO index 1 and defined by mode j_{0}.
    #                                           Since it is unique you do not have to perform a hermiticity or permutation relation.
    #                                           ...
    #                                           (Row 54) [4,20,0,0],[6,19,0,0],[7,18,0,0] = The entire row corresponds to a 3rd tier ADO with index 54, defined by 
    #                                           {j_{7},j_{6},j_{4}}. In the first element, [4,20,0,0], you remove mode j_{4} to connect to the 2nd tier ADO defined by 
    #                                           {j_{7},j_{6}}, with ADO index 20. In the second element, [6,19,0,0], you remove mode j_{4} to connect to the 2nd tier ADO defined by 
    #                                           {j_{7},j_{6}}, with ADO index 19. In the third element, [7,18,0,0], you remove mode j_{4} to connect to the 2nd tier ADO defined by 
    #                                           {j_{7},j_{6}}, with ADO index 18. Since these all connect to unique ADOs, there are no permutations or hermiticity relations required.
    

    def Gen_Ind_Minus(self):

        self.Index_minus = [[[]]]*(self.Hier_ind[self.Nmax,1]+1)                    # Initialize empty Index_minus array
        self.Index_minus[0] = [[-1]*4]*self.Nmax                                    # Fill the first row (zeroth tier) with -1

        ind_j = 1                                                                   # Start index of 1st tier ADOs
        for itrj1 in range(0,self.Nmodes,2):                                        # Loop through unique (even) 1st tier indices
            Index_minus_n = [[]]*self.Nmax                                          # Generate empty row of Index_minus
            Index_minus_n[0] = [itrj1,0,0,0]                                        # Indicates that this ADO connects to the 0th tier by removing mode j_{itrj1}
            for itrjrest in range(1,self.Nmax):                                     # Make the rest of the elements of this row -1
                Index_minus_n[itrjrest] = [-1]*4
            self.Index_minus[ind_j] = Index_minus_n                                 # Fill Index_minus with this row 
            ind_j += 1                                                              

        for itrn in range(2,self.Nmax+1):                                           # Loop through remaining tiers
            for itrjn in range(self.Hier_ind[itrn,0],self.Hier_ind[itrn,1]+1):      # Loop through all ADO indices of this tier
                j_old = self.Un_Ind[itrjn,0:itrn:1]                                 # Generate modes of this ADO from Un_Ind
                if itrn < self.Nmax:                                                # Check if we are below the maximum tier
                    Index_minus_n = [[]]*(itrn+1)                                   # If yes, then generate empty row
                    Index_minus_n[itrn:self.Nmax] = [[-1]*4]*(self.Nmax-itrn)       # and fill the remaining columns greater than itrn with -1
                else:                                       
                    Index_minus_n = [[]*4]*self.Nmax                                # If this is the maximum tier, generate empty row for all columns
                for itrjnm1 in range(itrn):                                         # Loop through all modes in ADO defined by j_old
                    j = j_old[itrjnm1]                                              # Identify mode to be removed
                    j_new = np.delete(j_old,itrjnm1)                                # Construct new ADO of tier itrn-1 with this mode removed
                    if j_new[0] %2 == 0:                                            # If the first mode in this new ADO is even, then it is one of the unique ADOs we explicitly include
                        minus_ind = 0                                               # Start index of new ADO at 0
                        for itrjnew in range(itrn-1):                               # Using the same process as in Index_plus, find the index of the new ADO by searching systematically in
                                                                                    # in Index_plus
                            minus_ind = self.Index_plus[minus_ind][int(j_new[itrjnew])][0]
                        minus_list = [int(j),minus_ind,0,0]                         # Generate element of Index_minus with the mode removed (int(j)) in the first position and the index of the new
                                                                                    # new ADO in the second position. No hermiticity relations or permutations are required, so other elements are 0
                        Index_minus_n[itrjnm1] = minus_list                         # Fill this row with the element in the appropriate place
                    elif j_new[0] %2 != 0:                                          # If the first mode in this new ADO is odd, however, then it is not one of the unique ADOs we 
                                                                                    # explicitly include, and we need to rewrite it using a hermiticity relation
                        ind_evodd = np.arange(itrn-1,dtype=int)                     # Similarly to the process in Gen_Ind_Plus, generate indices of modes in new ADO
                        ind_ev = ind_evodd[j_new%2==0]                              # Find indices of even modes
                        ind_odd = ind_evodd[j_new%2!=0]                             # Find indices of odd modes
                        j_new[ind_ev] += 1; j_new[ind_odd] -= 1                     # Perform conjugation
                        nperm = 0                                                   # Start permutation count at 0
                        if itrn > 2:                                                # If itrn = 2 then removing 1 mode leaves a 1st tier ADO with 1 mode, so no permutations are necessary
                            for itrpair in range(0,itrn-2):                         # Loop through pairs of modes in new ADO
                                if j_new[itrpair+1] < j_new[itrpair]:               # If pair is out of order, switch them and add 1 to number of permutations nperm
                                    j_new[[itrpair,itrpair+1]] = j_new[[itrpair+1,itrpair]]
                                    nperm += 1
                        minus_ind = 0                                               # Use same process to quickly find index of new ADO 
                        for itrjnew in range(itrn-1):
                            minus_ind = self.Index_plus[minus_ind][int(j_new[itrjnew])][0] 
                        minus_list = [int(j),minus_ind,nperm,1]                     # Create element of this row with mode to be removed, index of new ADO, number of permutations, and 1 to indicate
                                                                                    # a hermiticity relation was used.
                        Index_minus_n[itrjnm1] = minus_list
                self.Index_minus[itrjn] = Index_minus_n                             # Fill this row of Index_minus with information

    def Imp_Crit(self,n,indjn,eta_vec,gamma_vec,max_V_km):
        gamma_vec = np.sqrt(-gamma_vec)
        factor1 = 1.0
        factor2 = 1.0
        for itrjn in range(n):
            j = self.Un_Ind[indjn,itrjn]
            K_j = self.KsigLm[j,0]
            sig_j = self.KsigLm[j,1]
            l_j = self.KsigLm[j,2]
            m_j = self.KsigLm[j,3]
            eta_j = eta_vec[K_j][m_j][l_j]
            gamma_j = gamma_vec[K_j][l_j][sig_j].real
            factor1 = factor1*(eta_j/gamma_j)
            if itrjn < (n-1):
                factor2_denom = 0.0
                for itrjmjn in range(itrjn+1):
                    jm = self.Un_Ind[indjn,itrjmjn]
                    K_jm = self.KsigLm[jm,0]
                    sig_jm = self.KsigLm[jm,1]
                    l_jm = self.KsigLm[jm,2]
                    gamma_jm = np.real(gamma_vec[K_jm,l_jm,sig_jm])
                    factor2_denom += gamma_jm
                factor2 = factor2*((max_V_km**2)/factor2_denom)
        return abs(factor1*factor2)
        # ,factor1,factor2
        
    def Check_Threshold(self,tol,eta_vec,gamma_vec,max_V_km):
        self.Hier_ind_filtered = np.empty((self.Nmax+1,2),int)
        self.Hier_ind_filtered[0:2,0:2] = self.Hier_ind[0:2,0:2]
        self.Un_Ind_filtered = np.full((self.Hier_ind[self.Nmax,1]+1,self.Nmax),-1,dtype=int)

        retain = np.full(self.Hier_ind[self.Nmax,1]+1,True)
        filtering = np.empty(self.Hier_ind[self.Nmax,1]+1,int)
        ind_filt = 1
        self.Hier_ind_filtered[0,0:1] = 0
        for itrn in range(1,self.Nmax+1):
            for indjn in range(self.Hier_ind[itrn,0],self.Hier_ind[itrn,1]+1):
                if itrn < 2:
                    retain[indjn] = True
                    filtering[indjn] = ind_filt
                    self.Un_Ind_filtered[ind_filt] = copy.deepcopy(self.Un_Ind[indjn])
                    ind_filt += 1
                else:
                    IC = self.Imp_Crit(itrn,indjn,eta_vec,gamma_vec,max_V_km)
                    if IC < tol:
                        retain[indjn] = False
                        filtering[indjn] = -1
                    else:
                        retain[indjn] = True
                        filtering[indjn] = ind_filt
                        self.Un_Ind_filtered[ind_filt,:] = self.Un_Ind[indjn,:]
                        ind_filt += 1
            self.Hier_ind_filtered[itrn,0] = self.Hier_ind_filtered[itrn-1,1]+1
            self.Hier_ind_filtered[itrn,1] = ind_filt

        self.Un_Ind_filtered = self.Un_Ind_filtered[0:ind_filt,:]

        self.Index_minus_filtered = [[[]]]*(self.Hier_ind_filtered[self.Nmax,1])
        self.Index_minus_filtered[0] = self.Index_minus[0]
        self.Index_plus_filtered = [[[]]]*(self.Hier_ind_filtered[self.Nmax-1,1])
        self.Index_plus_filtered[0] = self.Index_plus[0]
        for itrn in range(1,self.Nmax+1):
            for indjn in range(self.Hier_ind[itrn,0],self.Hier_ind[itrn,1]+1):
                if retain[indjn] == True:
                    index_minus_temp = copy.deepcopy(self.Index_minus[indjn])
                    for itrtier in range(itrn):
                        ind_minus = index_minus_temp[itrtier][1]
                        if (ind_minus != -1): 
                            index_minus_temp[itrtier][1] = filtering[ind_minus]
                    self.Index_minus_filtered[filtering[indjn]] = index_minus_temp
                    if itrn < self.Nmax:
                        index_plus_temp = copy.deepcopy(self.Index_plus[indjn])
                        for itrmodes in range(self.Nmodes):
                            ind_plus = index_plus_temp[itrmodes][0]
                            if (ind_plus != -1): 
                                index_plus_temp[itrmodes][0] = filtering[ind_plus]
                        self.Index_plus_filtered[filtering[indjn]] = index_plus_temp
        return filtering,retain

    # ---------------------------------------------------------------------
    #
    #                   RETURN INFORMATION ON INDICES
    #
    # ---------------------------------------------------------------------
    #
    # This function returns the unfiltered index information calculated above.
    #
    # USAGE: KsigLm,Un_Ind,Hier_Ind,Index_minus,Index_plus,len_un_ind,len_index_plus,tier_index = Indices.Print_Ind_Info()
    

    def Print_Ind_Info(self):
        if self.wbl_YN == 0:
            return self.KsigLm,self.Un_Ind,self.Hier_ind,self.Index_minus,self.Index_plus,self.len_un_ind,self.len_index_plus,self.tier_index
        elif self.wbl_YN == 1:
            return self.KsigLm,self.Ksig0m,self.Un_Ind,self.Hier_ind,self.Index_minus,self.Index_plus,self.len_un_ind,self.len_index_plus,self.tier_index

    # ---------------------------------------------------------------------
    #
    #                RETURN INFORMATION ON FILTERED INDICES
    #
    # ---------------------------------------------------------------------
    #
    # This function returns the filtered index information calculated above.
    #
    # USAGE: KsigLm,Un_Ind_filtered,Hier_Ind_filtered,Index_minus_filtered,Index_plus_filtered = Indices.Print_Filtered_Ind_Info()
    

    def Print_Filtered_Ind_Info(self,tol,eta_vec,gamma_vec,max_V_km):
        self.Check_Threshold(tol,eta_vec,gamma_vec,max_V_km)
        if self.wbl_YN == 0:
            return self.KsigLm,self.Un_Ind_filtered,self.Hier_ind_filtered,self.Index_minus_filtered,self.Index_plus_filtered
        elif self.wbl_YN == 1:
            return self.KsigLm,self.Ksig0m,self.Un_Ind_filtered,self.Hier_ind_filtered,self.Index_minus_filtered,self.Index_plus_filtered


    # ---------------------------------------------------------------------
    #
    #                FIND LOCATION OF ADO GIVEN ITS MODES
    #
    # ---------------------------------------------------------------------
    #
    # USAGE: ind_j = Indices.Ind_Loc(j)
    #
    # INPUT:
    #       j                               Array containining modes defining ADO you want to know the index of
    #
    # OUTPUT:
    #       ind_j                           Corresponding ADO index

    def Ind_Loc(self,j):
        ind_j = 0
        for itrj in range(len(j)):
            ind_j = self.Index_plus[ind_j][int(j[itrj])][0]
        return ind_j

    # ---------------------------------------------------------------------
    #
    #          FIND INDEX OF NEW ADO GIVEN OLD ADO AND A NEW MODE
    #
    # ---------------------------------------------------------------------
    #
    # USAGE: ind_jplus = Indices.Ind_Loc_Plus(ind_j_old,j)
    #
    # INPUT:
    #        ind_j_old                      Index of old ADO
    # 
    #        j                              New mode to be added on to define new ADO
    #
    # OUTPUT:
    #       ind_jplus                       Corresponding index of newly created ADO

    def Ind_Loc_Plus(self,ind_j_old,j):
        ind_jplus = self.Index_plus[ind_j_old][int(j)][0]
        return ind_jplus

    # ---------------------------------------------------------------------
    #
    #        FIND INDEX OF NEW ADO GIVEN OLD ADO AND MODE TO REMOVE
    #
    # ---------------------------------------------------------------------
    #
    # USAGE: ind_jminus = Indices.Ind_Loc_Minus(ind_j_old,jrem)
    #
    # INPUT:
    #        ind_j_old                      Index of old ADO
    # 
    #        ind_jrem                       Position of mode to be removed
    #
    # OUTPUT:
    #       ind_jminus                       Corresponding index of newly created ADO

    def Ind_Loc_Minus(self,ind_j_old,ind_jrem):
        ind_jminus = self.Index_minus[ind_j_old][int(ind_jrem)][1]
        return ind_jminus

if __name__ == "__main__":                                                          # Produce example if run from command line of terminal
    
    Lmax = 1                                                                        # Define parameters
    Nmax = 3
    Nel = 1
    Nleads = 2
    Nsign = 2
    Nmodes = (Lmax+1)*Nel*Nleads*Nsign

    Index = Hierarchy_index(Nmax,Nel,Lmax,Nleads,Nsign,Nmodes)                      # Create object of Hierarchy_index class named 'Index' with example parameters
    KsigLm,Un_Ind,Hier_ind,Index_minus,Index_plus,len_un_ind,len_index_plus,tier_index = Index.Print_Ind_Info()
                                                                                    # Return various index information from object
    # Hierarchy_filtered = Hierarchy.Check_Threshold()                              # Uncomment (comment) if you (do not) want to apply filtering
    
    output_KsigLm = open('KsigLm.txt','w')                                          # Write the KsigLm array into a text file to view the different modes
    for itrj in range(0,Nmodes):
                output_KsigLm.write(''+str(KsigLm[itrj,:])+'\n')
    output_KsigLm.close()

    output_unique_indices = open('Unique_Indices.txt','w')                          # Write the Un_Ind array into a text file to view the unique indices for each tier
    for itrn in range(0,Nmax+1):
        if itrn == 0:
            output_unique_indices.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Un_Ind[0,:])+'\n')
        else:
            for itrjn in range(Hier_ind[itrn,0],Hier_ind[itrn,1]+1):
                output_unique_indices.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Un_Ind[itrjn,:])+'\n')
    output_unique_indices.close()

    output_index_plus = open('Index_plus.txt','w')                                  # Write the Index_plus array into a text file to view the connections between n and n+1 tiers
    for itrn in range(0,Nmax):
        if itrn == 0:
            output_index_plus.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Index_plus[0])+'\n')
        else:
            for itrjn in range(Hier_ind[itrn,0],Hier_ind[itrn,1]+1):
                output_index_plus.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Index_plus[itrjn])+'\n')
    output_index_plus.close()
    output_index_minus = open('Index_minus.txt','w')                                # Write the Index_minus array into a text file to view the connections between n and n-1 tiers
    for itrn in range(0,Nmax+1):
        if itrn == 0:
            output_index_minus.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Index_minus[0])+'\n')
        else:
            for itrjn in range(Hier_ind[itrn,0],Hier_ind[itrn,1]+1):
                output_index_minus.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Index_minus[itrjn])+'\n')
    output_index_minus.close()

    ############### UNCOMMENT CODE BELOW IF YOU WANT TO SEE FILTERED INDICES INFORMATION ####################

    # output_unique_indices_filtered = open('Unique_Indices_filtered.txt','w')
    # for itrn in range(0,Nmax+1):
    #     if itrn == 0:
    #         output_unique_indices_filtered.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Un_Ind_filtered[0,:])+'\n')
    #     else:
    #         for itrjn in range(Hier_ind_filtered[itrn,0],Hier_ind_filtered[itrn,1]+1):
    #             output_unique_indices_filtered.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Un_Ind_filtered[itrjn,:])+'\n')
    # output_unique_indices_filtered.close()

    # output_index_plus_filtered = open('Index_plus.txt_filtered','w')
    # for itrn in range(0,Nmax):
    #     if itrn == 0:
    #         output_index_plus_filtered.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Index_plus_filtered[0])+'\n')
    #     else:
    #         for itrjn in range(Hier_ind_filtered[itrn,0],Hier_ind_filtered[itrn,1]+1):
    #             output_index_plus_filtered.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Index_plus_filtered[itrjn])+'\n')
    # output_index_plus_filtered.close()

    # output_index_minus_filtered = open('Index_minus_filtered.txt','w')
    # for itrn in range(0,Nmax+1):
    #     if itrn == 0:
    #         output_index_minus_filtered.write('nth_tier['+str(0)+']\t= '+str(0)+'\tIndex['+str(0)+']\t= '+str(Index_minus_filtered[0])+'\n')
    #     else:
    #         for itrjn in range(Hier_ind_filtered[itrn,0],Hier_ind_filtered[itrn,1]+1):
    #             output_index_minus_filtered.write('nth_tier['+str(itrjn)+']\t= '+str(itrn)+'\tIndex['+str(itrjn)+']\t= '+str(Index_minus_filtered[itrjn])+'\n')
    # output_index_minus_filtered.close()