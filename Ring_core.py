from scipy.constants import c, hbar, pi
from scipy.fft import ifft, fft
from scipy.integrate import solve_ivp
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

class Ring:
    def __init__(self, ring_parameters) -> None:
        self.N = ring_parameters['N']
        self.n0 = ring_parameters['n0']
        self.n2 = ring_parameters['n2']
        self.FSR = ring_parameters['FSR']
        self.lambda0 = ring_parameters['lambda0']
        self.kappa = ring_parameters['kappa']
        self.eta = ring_parameters['eta']
        self.Veff = ring_parameters['Veff']
        self.D2 = ring_parameters['D2']
        self.Pin = ring_parameters['Pin']

        # Calculate further parameters
        self.Tr = 1 / self.FSR # Roundtrip time [s]
        self.freq0 = c / self.lambda0 # CW pump frequency [Hz]
        self.omega0 = 2 * pi * self.freq0 # CW pump angular frequency [rad/s]
        self.g = hbar * self.omega0**2 * c * self.n2 / (self.n0**2 * self.Veff) # Nonlinear coupling coefficient
        self.kappa = 2 * pi * self.kappa # Optical linewidth [rad/s]
        self.f = np.zeros(self.N) # Normalized pump field vector
        self.f[int((self.N / 2) + 1)] = np.sqrt((8 * self.g * self.eta / self.kappa**2) * (self.Pin / (hbar * self.omega0))) # Normalized pump field
        self.D2 = 2 * pi * self.D2 # Second order dispersion [rad/s]
        self.d2 = (2 / self.kappa) * self.D2 # Normalized second order dispersion
        self.mu = np.arange(-(self.N - 1) / 2, ((self.N - 1) / 2) + 1, 1) # Mode numbers relative to the pumped mode
        self.dint = (self.d2 / 2) * self.mu**2 # Normalized integrated dispersion
    

    def num_sim(self, parameters, simulation_options):
        # Retrieve simulation options
        Effects = simulation_options['Effects']
        Noise = simulation_options['Noise']

        # Retrieve simulation parameters
        dseta_start = parameters['dseta_start']
        dseta_end = parameters['dseta_end'] 
        dseta_step = parameters['dseta_step']
        roundtrips_step = parameters['roundtrips_step']
        Amu0 = parameters['Amu0']
        if Effects == 'Thermal':
            theta0 = parameters['theta0']
            tauT = parameters['tauT']
            n2T = parameters['n2T']
        elif Effects == 'Avoided mode crossings':
            theta0, tauT, n2T = 0, 0.1, 0
            mode_perturbated = parameters['mode_perturbated']
        else:
            theta0, tauT, n2T = 0, 0.1, 0

        # Calculate further parameters
        dseta = np.arange(dseta_start, dseta_end + dseta_step, dseta_step) # Normalized detuning (vector)
        tau_step = (self.kappa / 2) * self.Tr * roundtrips_step # Normalized time per tuning step
        amu = np.zeros((len(dseta), self.N), dtype=np.complex_) # 2D array to store normalized fields
        amu[0, :] = ifft(np.sqrt(2 * self.g / self.kappa) * Amu0) # Store normalized initial field
        theta = np.zeros(len(dseta), dtype=np.complex_)
        theta[0] = theta0
        aux = 1j if Effects == 'Thermal' else 0
        if Effects == 'Avoided mode crossings':
            self.dint[((self.N + 1) / 2) + mode_perturbated] = 0

        # Iterate over all detuning values
        for i in tqdm(range(len(dseta) - 1)):
            dseta_current = dseta[i] # Current detuning value
            y0 = np.insert(amu[i, :], self.N, theta[i]) # Initial conditions to solve LLE
            if i == 0:
                def LLE(tau, y):
                    amu_ = y[:-1]
                    theta_ = y[-1]
                    damu_ = -(1 + 1j * (dseta_current + dseta_step * (tau / tau_step) + self.dint) - aux * theta_) * amu_ + 1j * ifft(np.abs(fft(amu_))**2 * fft(amu_)) + self.f # LLE
                    dtheta_ = (2 / (self.kappa * tauT)) * ((n2T / self.n2) * np.sum(np.abs(amu_)**2) - theta_)
                    dy = np.insert(damu_, self.N, dtheta_)
                    return dy 
            solution = solve_ivp(LLE, [0, tau_step], y0) # Solve LLE
            noise = ifft(np.sqrt(2 * self.g / self.kappa) * (np.random.randn(self.N) + np.random.randn(self.N) * 1j)) if Noise else 0
            amu[i + 1, :] = solution.y[:-1, -1] + noise # Store field solution
            theta[i + 1] = solution.y[-1, -1]
        
        return dseta, amu, theta


    def plot_res(self, dseta_forward, amu_forward, dseta_backward = np.array([]), amu_backward = np.array([]), dseta_snap=0):
        fig = plt.figure(figsize=(14, 14))
        plt.rcParams['font.size'] = '18'

        avg_intr_int_fwrd = np.sum(np.abs(amu_forward)**2, 1) # Average intracavity intensity of forward tuning
        ax1 = plt.subplot(311)
        ax1.plot(dseta_forward, avg_intr_int_fwrd, label='Forward tuning', linewidth=2)
        ax1.axvline(dseta_snap, color='black', linestyle='--', linewidth=2, label='$\zeta$ = {}'.format(dseta_snap)) # Plot dseta_snap
        ax1.set_xlabel('Normalized detuning, $\zeta$')
        ax1.set_ylabel('Average intracavity intensity [a. u.]')
        ax1.set_xlim(np.min(dseta_forward), np.max(dseta_forward))
        ax1.set_ylim(0, np.max(avg_intr_int_fwrd) * 1.05)
        ax1.legend()

        ax2 = plt.subplot(312) if amu_backward.size == 0 else plt.subplot(323)
        self.plot_optical_spectrum(dseta_forward, amu_forward, dseta_snap, ax2)

        ax4 = plt.subplot(313) if amu_backward.size == 0 else plt.subplot(325)
        self.plot_intracavity_power(dseta_forward, amu_forward, dseta_snap, ax4)

        if amu_backward.size != 0:
            ax3 = plt.subplot(324)
            self.plot_optical_spectrum(dseta_backward, amu_backward, dseta_snap, ax3)

            ax5 = plt.subplot(326)
            self.plot_intracavity_power(dseta_backward, amu_backward, dseta_snap, ax5)
            
        fig.canvas.header_visible = False
        plt.tight_layout()
        plt.show()


    def plot_optical_spectrum(self, dseta, amu, dseta_snap, ax):
        index = np.abs(dseta - dseta_snap).argmin()
        omegamu = self.omega0 + 2 * pi * self.FSR * self.mu + (2 * pi * self.D2 * self.mu**2 / 2)
        freqmu = omegamu / (2 * pi) # Resonance frequencies [Hz]
        fin = np.max(self.f) * np.ones(self.N)
        fout = self.f - 2 * self.eta * np.abs(amu[index, :])
        optical_spectrum = (np.abs(fout)**2 / np.abs(fin)**2) * self.Pin
        optical_spectrum = 10 * np.log10(1000 * optical_spectrum)

        ax.stem(freqmu * 1e-12, optical_spectrum, markerfmt=",", bottom=-30)
        ax.set_xlabel('Resonance frequencies [THz], $f_\mu$')
        ax.set_ylabel('Optical spectrum [dBm]')
        ax.set_xlim(freqmu[0] * 1e-12, freqmu[-1] * 1e-12)
        ax.set_ylim(-30, 10 * np.log10(1000 * self.Pin) * 1.05)


    def plot_intracavity_power(self, dseta, amu, dseta_snap, ax):
        index = np.abs(dseta - dseta_snap).argmin()
        ring_circumference = np.linspace(-pi, pi, 4096)
        temp = np.zeros_like(ring_circumference, dtype=np.complex_)
        temp[np.int32(self.mu + (4096 / 2) + 1)] = amu[index, :]
        intracavity_power = np.abs(fft(temp))**2

        ax.plot(ring_circumference, intracavity_power, linewidth=2)
        ax.set_xlabel('Ring circumference [rad], $\phi$')
        ax.set_ylabel('Intracavity power [a. u]')
        ax.set_xlim(-pi, pi)
        ax.set_ylim(0, np.max(intracavity_power) * 1.05)
        ax.set_xticks([-pi, -pi/2, 0, pi/2, pi])
        ax.set_xticklabels(['$-\pi$', '$-\pi$/2', '0', '$\pi$/2', '$\pi$'])
