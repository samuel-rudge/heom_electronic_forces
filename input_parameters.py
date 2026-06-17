import numpy as np
from constants import * # pylint disable=unused-wildcard-import
import const_ARK
import CreAnn
import scipy
import math

# Define system constraints

Nel = 2
N_qu_vib_modes = 0
N_cl_vib_modes = 1
max_occ_qu_vib_modes = [0]                                                                # Maximum number of phonons allowed in each mode
dim_vib_mode_qu = np.prod(np.array(max_occ_qu_vib_modes,dtype=int)+1)                             # Number of bosonic Fock states
dim_el = 2**Nel                                                                     # Number of fermionic Fock states
dim_rho = dim_el*dim_vib_mode_qu
N_el_vib_int_qu = 0
N_el_vib_int_cl = 1
small_polaron_yn = 0

# Define hierarchical constraints

Nmax = 2
Npoles_pade = 10
Npoles_barycentric = 15
Nsupport_points_barycentric = 2000
tol_F = 1e-5
tol_Gamma_barycentric = 1e-5
tol_fermi_symmetrized_barycentric = 1e-3
Nleads = 2 
Nsign = 2                                                                         # Types of second-quantization operators. Always = 2: annihilation (d) and creation (ddag).
filtering_YN = 0
if filtering_YN == 1:
    tol = 0.0

# Define constraints of bath spectral functions

wbl_YN = 0
analytic_spectral_function_decomposition = True
pole_choice = "barycentric"
if wbl_YN == 0:                                                                    # If not under the wide-band limit, they are assumed to have a Lorentzian density of states
    Nmodes_pade = (Npoles_pade+1)*Nel*Nleads*Nsign                                            # Calculate number of modes outside of wide-band limit
    Nmodes_barycentric = (Npoles_barycentric+1)*Nel*Nleads*Nsign                                            # Calculate number of modes outside of wide-band limit
    if Nleads == 2:
        specwidth = np.array([10.0,10.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
        symmetrized_fermi_specwidth = np.array([10.0,10.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
    if Nleads == 1:
        specwidth = np.array([10.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
        symmetrized_fermi_specwidth = np.array([10.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
    Ncutoff = 1
elif wbl_YN == 1:
    Nmodes_pade = Npoles_pade*Nel*Nleads*Nsign                                            # Calculate number of modes outside of wide-band limit
    Nmodes_barycentric = Npoles_barycentric*Nel*Nleads*Nsign                                            # Calculate number of modes outside of wide-band limit
    if Nleads == 2:
        specwidth = np.array([5.0,5.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
        symmetrized_fermi_specwidth = np.array([5.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
    if Nleads == 1:
        specwidth = np.array([10.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
        symmetrized_fermi_specwidth = np.array([5.0],dtype=float)                                 # Set spectral width of bath Lorentzians outside of wide-band limit
    Ncutoff = 1
    Nwblmodes = Nel*Nleads*Nsign

# Define propagation and derivative constraints

dt_init = 1e-1
dt_min = 1e-2
atol = 1e-6
rtol = 1e-6
fac = 0.38**(1/5)
facmin = 0.2
facmax_init = 10.0
rk_coeff,rk_coeffhat = const_ARK.dvecs('Dormand Prince')
max_expan_order = np.shape(rk_coeff)[1]                                         # Maximum expansion order for Euler's solution of DE
tol_sssd = 1e-2
tol_heom_prop = 1e-12
tol_sssd_prop = 1e-12
tol_friction = 1e-12
tol_corrfunc = 1e-12
# rho_0 = np.diag(np.ones(dim_rho,dtype=complex))/dim_rho
rho_0 = np.zeros((dim_rho,dim_rho),dtype=complex) ; rho_0[0,0] = 1.
dim_krylov_space = 15000
tol_gmres_ss = 1e-10
nthreads_liouvillian = 1
nthreads_x_grid = 19
printing_timestep = 0.1
checking_partitions_and_voltages = False
parallelize_x_grid = False
checking_timesteps_yn = False
redo_everything = False
print_integrand_yn = False
checking_time_interval = 300
tol_checking_time = 1e-6
n_prony_terms = 500

max_time = 1500
min_time = 500
# Define bath parameters

voltage = 0.4
dv = 0.05
if Nleads == 2:
    muvec = np.array([voltage/2,-voltage/2],dtype=float)                      # Chemical potentials of left (first) and right (second) electrodes
elif Nleads == 1:
    muvec = np.array([voltage/2],dtype=float)                      # Chemical potentials of left (first) and right (second) electrodes

Kelvin_T = 300
Temp = Kelvin_T*k_B                                                      # Electrode temperature, assumed to be same for both electrodes

# Define system classical vibrational parameters

freq_vector_cl_vib_modes = [0.03]
el_vib_int_cl = [0.01]
Vib_Freq_cl = np.zeros(N_cl_vib_modes,dtype=float)
El_Nuclear_Couplings_cl = np.zeros(N_cl_vib_modes,dtype=float)
small_polaron_shift_cl = np.zeros(N_cl_vib_modes,dtype=float)
for itr_cl_vib_modes in range(N_cl_vib_modes):
    Vib_Freq_cl[itr_cl_vib_modes] = freq_vector_cl_vib_modes[itr_cl_vib_modes]
    El_Nuclear_Couplings_cl[itr_cl_vib_modes] = np.sqrt(2)*el_vib_int_cl[itr_cl_vib_modes]
    small_polaron_shift_cl[itr_cl_vib_modes] = (el_vib_int_cl[itr_cl_vib_modes]**2)/Vib_Freq_cl[itr_cl_vib_modes]

# Define system quantum vibrational parameters

freq_vector_qu_vib_modes = [0.3]
el_vib_int_qu = [0.45]
Vib_Freq_qu = np.zeros(N_qu_vib_modes,dtype=float)
El_Nuclear_Couplings_qu = np.zeros(N_qu_vib_modes,dtype=float)
small_polaron_shift_qu = np.zeros(N_qu_vib_modes,dtype=float)
for itr_qu_vib_modes in range(N_qu_vib_modes):
    Vib_Freq_qu[itr_qu_vib_modes] = freq_vector_qu_vib_modes[itr_qu_vib_modes]
    El_Nuclear_Couplings_qu[itr_qu_vib_modes] = el_vib_int_qu[itr_qu_vib_modes]
    small_polaron_shift_qu[itr_qu_vib_modes] = (El_Nuclear_Couplings_qu[itr_qu_vib_modes]**2)/Vib_Freq_qu[itr_qu_vib_modes]

# Define system electronic parameters

delta_parameter = 0.2
energies_vector = [delta_parameter/2,-delta_parameter/2]
hopping_vector = [0.0,0.0]
elel_interaction_vector = [0.0]
degenerate_levels = True
Single_El_Int = np.zeros((Nel,Nel),dtype=float)        # Energies of levels included in system, as well as hopping between levels
Double_El_Int = np.zeros((Nel,Nel),dtype=float)        # Coulomb interactions between fermions - expressed as lower triangular matrix filled with two-particle interactions
hopping_count = 0
if Nel == 1:
    Single_El_Int[0,0] = energies_vector[0] #+ small_polaron_shift_cl[0]
    #    Single_El_Int[0,0] = energies_vector[0] + small_polaron_shift_qu[0]
elif Nel > 1:
    for itr_el1 in range(Nel):
        Single_El_Int[itr_el1,itr_el1] = energies_vector[itr_el1] #+ small_polaron_shift_cl[0]
        for itr_el2 in range(itr_el1+1,Nel):
            Single_El_Int[itr_el1,itr_el2] = hopping_vector[hopping_count]
            Single_El_Int[itr_el2,itr_el1] = hopping_vector[hopping_count]
            Double_El_Int[itr_el1,itr_el2] = elel_interaction_vector[hopping_count] #+ small_polaron_shift_cl[0]
            Double_El_Int[itr_el2,itr_el1] = elel_interaction_vector[hopping_count] #+ small_polaron_shift_cl[0]
            hopping_count+=1

# Define vibrational vector 

max_energy_diff = 0.45
dx_grid = 1
x_min_total = 0.001
x_max_total = x_min_total+dx_grid*1.001
len_x_vec_total = int((x_max_total - x_min_total)/dx_grid)

if (len_x_vec_total > 1):
    dx_init = (x_max_total - x_min_total)/len_x_vec_total
else:
    dx_init = 1.

n_x_partitions = 1
len_x_vec = math.ceil(len_x_vec_total/n_x_partitions)
x_partition_count = 0
ind_x_total = np.arange(x_partition_count*len_x_vec,(x_partition_count+1)*len_x_vec)
np.save("ind_x_total.npy",ind_x_total)

x_min = x_min_total + len_x_vec*x_partition_count*dx_grid
x_max = x_min_total + len_x_vec*(x_partition_count+1)*dx_grid-dx_grid
x_vec = np.linspace(x_min,x_max,len_x_vec)

dx_sssd = 0.001
rho_0_x_arr = np.repeat(rho_0[:, :, np.newaxis], len_x_vec, axis=2)
ind_x_0 = np.argmin(np.abs(x_vec))
U_0_vec = (0.5*Vib_Freq_cl[0]*x_vec**2)
U_1_vec = U_0_vec + (Single_El_Int[0,0] + El_Nuclear_Couplings_cl[0]*x_vec)
dx_vec = np.array([0,2*dx_sssd,dx_sssd,-dx_sssd,-2*dx_sssd],dtype=float)

# DEFINE GAMMA PARAMETERS ###

Gamma_choice = 25/1000
V_Km = np.sqrt(Gamma_choice/(2*np.pi))
# el_lead_couplings = np.full((Nleads,Nel,len_x_vec,5),V_Km,dtype=float)
el_lead_couplings = np.zeros((Nleads,Nel,len_x_vec,5),dtype=float)
dV_dx = np.zeros((Nleads,Nel,len_x_vec),dtype=float)

for itrleads in range(Nleads):
    for itrel in range(Nel):
        for itrx in range(len_x_vec):
            for itrdx in range(5):
                el_lead_couplings[0,0,itrx,itrdx] = V_Km
                el_lead_couplings[1,1,itrx,itrdx] = V_Km

