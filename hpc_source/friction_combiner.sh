#!/bin/bash
code_dir="omega_qu_300meV_lambda_qu_450meV_final_v6"
friction_file="friction.dat"
friction_mol_file="friction_mol.dat"
friction_molleads_file="friction_molleads.dat"
corrfunc_file="corrfunc.dat"
corrfunc_mol_file="corrfunc_mol.dat"
corrfunc_molleads_file="corrfunc_molleads.dat"
adiabatic_force_file="adiabatic_force.dat"
adiabatic_force_mol_file="adiabatic_force_mol.dat"
adiabatic_force_molleads_file="adiabatic_force_molleads.dat"
current_file="current.dat"
average_vib_excitation_diag_file="average_vib_excitation_diag.dat"
average_vib_excitation_nondiag_file="average_vib_exctiation_nondiag.dat"

for friction_dir in $code_dir/*/*/*/*/; do
    cd $friction_dir
    if test -f $friction_file; then
        rm $friction_file
	rm $friction_mol_file
	rm $friction_molleads_file
    fi
    if test -f $corrfunc_file; then
        rm $corrfunc_file
	rm $corrfunc_mol_file
	rm $corrfunc_molleads_file
    fi
    if test -f $adiabatic_force_file; then
        rm $adiabatic_force_file
	rm $adiabatic_force_mol_file
	rm $adiabatic_force_molleads_file
    fi
    if test -f $current_file; then
        rm $current_file
    fi
    if test -f $average_vib_excitation_diag_file; then
        rm $average_vib_excitation_diag_file
	rm $average_vib_excitation_nondiag_file
    fi
    touch $friction_file
    echo -e "Coord. Val.\tFriction" > $friction_file
    touch $friction_mol_file
    echo -e "Coord. Val.\tFriction" > $friction_mol_file
    touch $friction_molleads_file
    echo -e "Coord. Val.\tFriction" > $friction_molleads_file
    touch $corrfunc_file
    echo -e "Coord. Val.\tCorr. Func." > $corrfunc_file
    touch $corrfunc_mol_file
    echo -e "Coord. Val.\tCorr. Func." > $corrfunc_mol_file
    touch $corrfunc_molleads_file
    echo -e "Coord. Val.\tCorr. Func." > $corrfunc_molleads_file
    touch $adiabatic_force_file
    echo -e "Coord. Val.\tAdiabatic Force" > $adiabatic_force_file
    touch $adiabatic_force_mol_file
    echo -e "Coord. Val.\tAdiabatic Force" > $adiabatic_force_mol_file
    touch $adiabatic_force_molleads_file
    echo -e "Coord. Val.\tAdiabatic Force" > $adiabatic_force_molleads_file
    touch $current_file
    echo -e "Coord. Val.\tLead 0\tLead 1" > $current_file
    touch $average_vib_excitation_diag_file
    echo -e "Coord. Val.\tVib. Exc." > $average_vib_excitation_diag_file
    touch $average_vib_excitation_nondiag_file
    echo -e "Coord. Val.\tVib. Exc." > $average_vib_excitation_nondiag_file


    for x_partition_dir in */ ; do
        echo $x_partition_dir
	cat $adiabatic_force_file <(tail -n +2 "$x_partition_dir/adiabatic_force.dat") >> $adiabatic_force_file
	cat $adiabatic_force_mol_file <(tail -n +2 "$x_partition_dir/adiabatic_force_mol.dat") >> $adiabatic_force_mol_file
	cat $adiabatic_force_molleads_file <(tail -n +2 "$x_partition_dir/adiabatic_force_molleads.dat") >> $adiabatic_force_molleads_file
	cat $friction_file <(tail -n +2 "$x_partition_dir/friction.dat") >> $friction_file
	cat $friction_mol_file <(tail -n +2 "$x_partition_dir/friction_mol.dat") >> $friction_mol_file
	cat $friction_molleads_file <(tail -n +2 "$x_partition_dir/friction_molleads.dat") >> $friction_molleads_file
	cat $corrfunc_file <(tail -n +2 "$x_partition_dir/corrfunc.dat") >> $corrfunc_file
	cat $corrfunc_mol_file <(tail -n +2 "$x_partition_dir/corrfunc_mol.dat") >> $corrfunc_mol_file
	cat $corrfunc_molleads_file <(tail -n +2 "$x_partition_dir/corrfunc_molleads.dat") >> $corrfunc_molleads_file
	cat $current_file <(tail -n +2 "$x_partition_dir/current.dat") >> $current_file
	cat $average_vib_excitation_diag_file <(tail -n +2 "$x_partition_dir/average_vib_excitation_diag.dat") >> $average_vib_excitation_diag_file
	cat $average_vib_excitation_nondiag_file <(tail -n +2 "$x_partition_dir/average_vib_excitation_nondiag.dat") >> $average_vib_excitation_nondiag_file
    done
    pwd
    cd -
done

