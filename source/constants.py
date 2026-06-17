from scipy.constants import physical_constants as phy_const

hartree = phy_const["Hartree energy in eV"][0]                            # 1 hartree=27.21138602 eV
hartree_to_K = phy_const["hartree-kelvin relationship"][0]                # 1 hartree=315775.13 K
hartree_to_cm_1 = phy_const["hartree-inverse meter relationship"][0]/100  # 1 hartree=219474.6313702 cm^{-1}
k_B = phy_const["Boltzmann constant in eV/K"][0]                          # k_B in units of eV/K
au_to_fs = phy_const["atomic unit of time"][0]*1e15                       # 1 atomic unit of time = 2.41888432e-2 fs
h_planck=phy_const["Planck constant in eV s"][0]
au_to_amp=phy_const["atomic unit of current"][0]

if __name__=='__main__':
    print('constants:\t k_B = '+str(k_B))
    print('constants:\t h_planck = '+str(h_planck))
    print('constants:\t 1/h_planck = '+str(1.0/(1.e-15/h_planck)))
    
