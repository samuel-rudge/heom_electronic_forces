import numpy as np
from matplotlib.lines import Line2D
import matplotlib.transforms as mtransforms
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import sys,glob,os
import matplotlib as mpl
mpl.rcParams['savefig.directory'] = os.getcwd()
from input_parameters import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes', linewidth=2)
plt.rc('text.latex', preamble=r'\boldmath')

dt_sample = dt_min/10

# Import and sample friction

time_vec_friction = np.genfromtxt("friction_integrand_heom.dat")[:,0]
time_vec_friction_uniform = np.arange(0,np.max(time_vec_friction),dt_sample)
N_timesteps_friction = len(time_vec_friction_uniform)
friction_integrand = -np.genfromtxt("friction_integrand_heom.dat")[:,1]
friction_integrand_interp = interp1d(time_vec_friction,friction_integrand)
friction_integrand_uniform = friction_integrand_interp(time_vec_friction_uniform)

friction_integrand_fft = np.fft.fft(friction_integrand_uniform)           # complex Fourier transform
freqs_friction = np.fft.fftfreq(N_timesteps_friction,dt_sample)
friction_integrand_fft_shifted = np.fft.fftshift(friction_integrand_fft)
freqs_friction_shifted = np.fft.fftshift(freqs_friction)


# Import and sample correlation function

time_vec_corrfunc = np.genfromtxt("corrfunc_integrand_heom.dat")[:,0]
time_vec_corrfunc_uniform = np.arange(0,np.max(time_vec_corrfunc),dt_sample)
N_timesteps_corrfunc = len(time_vec_corrfunc_uniform)
corrfunc_integrand = np.genfromtxt("corrfunc_integrand_heom.dat")[:,1]
corrfunc_integrand_interp = interp1d(time_vec_corrfunc,corrfunc_integrand)
corrfunc_integrand_uniform = corrfunc_integrand_interp(time_vec_corrfunc_uniform)

corrfunc_integrand_fft = np.fft.fft(corrfunc_integrand_uniform)           # complex Fourier transform
freqs_corrfunc = np.fft.fftfreq(N_timesteps_corrfunc,dt_sample)
corrfunc_integrand_fft_shifted = np.fft.fftshift(corrfunc_integrand_fft)
freqs_corrfunc_shifted = np.fft.fftshift(freqs_corrfunc)

# Plotting

fig, ax = plt.subplots()
ax.plot(freqs_friction_shifted*2*np.pi,np.real(friction_integrand_fft_shifted*dt_sample),
        color="b",linewidth=2,label=r"$\displaystyle Re\tilde{\gamma}(E)$")
# ax.plot(freqs_corrfunc_shifted*2*np.pi,np.real(corrfunc_integrand_fft_shifted*dt_sample),
#         color="r",linewidth=2,label=r"$\displaystyle Re\tilde{D}(E)$")
ax.plot(freqs_friction_shifted*2*np.pi,np.imag(friction_integrand_fft_shifted*dt_sample),
        color="b",linestyle="--",linewidth=2,label=r"$\displaystyle Im\tilde{\gamma}(E)$")
# ax.plot(freqs_corrfunc_shifted*2*np.pi,np.imag(corrfunc_integrand_fft_shifted*dt_sample),
#         color="r",linestyle="--",linewidth=2,label=r"$\displaystyle Im\tilde{D}(E)$")
ax.set_xlabel(r"$\displaystyle E \: [\mbox{\textbf{eV}}] $",color='black',fontsize=24,fontweight='bold')
ax.set_ylabel(r"$\displaystyle \mbox{\textbf{Fourier Transform}}$",color='black',fontsize=24,fontweight='bold')
ax.tick_params(axis='x', labelcolor='black',length=6, width=2,labelsize=20)
ax.tick_params(axis='y', labelcolor='black',length=6, width=2,labelsize=20)
handles,labels = ax.get_legend_handles_labels()
voltage_label = ax.legend(handles,labels,loc='upper right',fontsize=16)
ax.set_xlim([-0.5*2*np.pi,0.5*2*np.pi])
# ax.set_ylim([0,np.max(np.array(spin_selectivity_arr))*1.05])
plt.tight_layout()
#plt.savefig(f"fourier_transform_correlation_and_friction_function_gamma_{Gamma_choice}eV_delta_{delta_parameter}eV_x_{x_min_total}.pdf")
plt.show()
