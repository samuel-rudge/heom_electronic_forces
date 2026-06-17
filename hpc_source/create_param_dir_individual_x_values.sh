#!/bin/bash
code_dir="2l1m_negative_friction_w_time_dependence"
if [ ! -d "$code_dir" ]
then
        mkdir $code_dir
fi

########################## BATH AND HEOM COMPUTATIONAL CONSTRAINT PARAMETERS ##########################

gamma_list=(100)
temp_list=(300)
voltage_list=(0.0 0.40)
ncores=20
tier_list=(2)
#3)

########################## TIME PROPAGATION PARAMETERS ###################################

dt_init="1e-1"
dt_min="1e-2"
wbl_YN=0
dim_krylov_space=15000
tol_friction="1e-12"
tol_corrfunc="1e-12"
max_time=10000
min_time=500
atol="1e-6"
rtol="1e-6"

########################## SYSTEM PARAMETERS ############################

x_min=-25
x_max=25
dx=0.1

x_list=(0 5 10)
max_energy_diff=0.40
nleads=2
nel=2
n_cl_vib_modes=1
n_cl_el_vib_int=1
delta_parameter_list=(0.2 0 -0.2)
#energies_vector="[0.1,-0.1]"
degenerate_levels="True"
hopping_vector="[0.0,0.0]"
elel_interaction_vector="[0.0]"
freq_vector_cl_vib_modes="[0.03]"
el_vib_int_cl="[0.01]"

working_directory="/home/sr1160/Documents/Postdoc/Projects/Project_negative_friction/Code/semiclassical_methods"
moving_files_directory="moving_files_individual_x_values/"

for itrgamma_temp in ${!gamma_list[@]}
 do
        gamma_value=${gamma_list[itrgamma_temp]}
        temp_value=${temp_list[itrgamma_temp]}
        gammatemp_dir="gamma_${gamma_value}meV_temp_${temp_value}K"
        if [ ! -d "$code_dir/$gammatemp_dir" ]
        then
                mkdir $code_dir/$gammatemp_dir
        fi
for itrdelta in ${!delta_parameter_list[@]}
 do
        delta_value=${delta_parameter_list[itrdelta]}
        delta_dir="delta_${delta_value}eV"
        if [ ! -d "$code_dir/$gammatemp_dir/$delta_dir" ]
        then
                mkdir $code_dir/$gammatemp_dir/$delta_dir
        fi
for itr_tier in ${!tier_list[@]}
 do
        tier_value=${tier_list[itr_tier]}
        tier_dir="tier_${tier_value}"
        if [ ! -d "$code_dir/$gammatemp_dir/$delta_dir/$tier_dir" ]
        then
                mkdir $code_dir/$gammatemp_dir/$delta_dir/$tier_dir
        fi	
    for itrvoltage in ${!voltage_list[@]}
     do
        voltage_value=${voltage_list[itrvoltage]}
        voltage_dir="voltage_${voltage_value}eV"
        if [ ! -d "$code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir" ]
        then
              mkdir $code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir
        fi
	x_value=$x_min
	while (( $(echo "$x_value <= $x_max" | bc -l) )); do
	        x_dir="x_${x_value}"
	        if [ ! -d "$code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir/$x_dir" ]
	        then
	              mkdir $code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir/$x_dir
		fi

        cp -r $working_directory/$moving_files_directory/. $code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir/$x_dir/
        cd $code_dir/$gammatemp_dir/$delta_dir/$tier_dir/$voltage_dir/$x_dir

	sed -i  "s%.*delta_parameter = .*%delta_parameter = $delta_value%" input_parameters.py
#        sed -i  "s%.*energies_vector = .*%energies_vector = $energies_vector%" input_parameters.py
        sed -i  "s%.*Nmax = .*%Nmax = $tier_value%" input_parameters.py
        sed -i  "s%.*freq_vector_cl_vib_modes = .*%freq_vector_cl_vib_modes = $freq_vector_cl_vib_modes%" input_parameters.py
        sed -i  "s%.*el_vib_int_cl = \[.*%el_vib_int_cl = $el_vib_int_cl%" input_parameters.py
#        sed -i  "s%.*dx_grid = .*%dx_grid = $dx_grid%" input_parameters.py
	sed -i  "s%.*x_min_total = .*%x_min_total = $x_value%" input_parameters.py
        sed -i  "s%.*max_energy_diff = .*%max_energy_diff = $max_energy_diff%" input_parameters.py
        sed -i  "s%.*Kelvin_T = .*%Kelvin_T = $temp_value%" input_parameters.py
        sed -i  "s%.*voltage = .*%voltage = $voltage_value%" input_parameters.py
        sed -i  "s%.*Gamma_choice = .*%Gamma_choice = $gamma_value/1000%" input_parameters.py
        sed -i  "s%.*dim_krylov_space = .*%dim_krylov_space = $dim_krylov_space%" input_parameters.py
        sed -i  "s%.*tol_friction = .*%tol_friction = $tol_friction%" input_parameters.py
        sed -i  "s%.*tol_corrfunc = .*%tol_corrfunc = $tol_corrfunc%" input_parameters.py
        sed -i  "s%.*max_time = .*%max_time = $max_time%" input_parameters.py
        sed -i  "s%.*min_time = .*%min_time = $min_time%" input_parameters.py
        sed -i  "s%.*nthreads_x_grid = .*%nthreads_x_grid = 1%" input_parameters.py
        sed -i  "s%.*nthreads_liouvillian = .*%nthreads_liouvillian = $ncores%" input_parameters.py
        python3 friction_heom_main.py ss
        sed -i  "s%.*nthreads_x_grid = .*%nthreads_x_grid = $ncores%" input_parameters.py
        sed -i  "s%.*nthreads_liouvillian = .*%nthreads_liouvillian = 1%" input_parameters.py
	x_value=$(echo "$x_value + $dx" | bc -l)
        python3 friction_heom_main.py markovian
#        python3 fourier_transform_and_plot.py
	cd -
      done
   done
done
done
done
