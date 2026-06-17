#!/bin/bash
code_dir="2l1m_negative_friction_w_time_dependence"

find $code_dir -maxdepth 6 -name 'friction.dat' -o -name 'friction_integrand_heom.dat' -o -name 'corrfunc.dat' -o -name 'corrfunc_integrand_heom.dat' -o -name 'current_ad.dat' -o -name 'current_na.dat' -o -name 'adiabatic_force.dat' | zip $code_dir.zip -@

