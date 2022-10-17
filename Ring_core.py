from scipy.constants import c, hbar, pi
from scipy.fft import ifft, fft
from scipy.integrate import solve_ivp
import numpy as np
from tqdm import tqdm

class Ring:
    def __init__(self, ring_param) -> None:
        self.N = ring_param['N']
        self.n0 = ring_param['n0']
        self.n2 = ring_param['n2']
        self.FSR = ring_param['FSR']
        self.lambda0 = ring_param['lambda0']
        self.kappa = ring_param['kappa']
        self.eta = ring_param['eta']
        self.Veff = ring_param['Veff']
        self.D2 = ring_param['D2']
        self.n2T = ring_param['n2T'] 
        self.Pin = ring_param['Pin']

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
    
    def num_sim(self, params, sim_opts):
        # Retrieve simulation parameters
        dseta_start = params['dseta_start']
        dseta_end = params['dseta_end'] 
        dseta_step = params['dseta_step']
        roundtrips_step = params['roundtrips_step']
        Amu0 = params['Amu0'] 

        # Calculate further parameters
        dseta = np.arange(dseta_start, dseta_end + dseta_step, dseta_step) # Normalized detuning (vector)
        tau_step = (self.kappa / 2) * self.Tr * roundtrips_step # Normalized time per tuning step
        amu = np.zeros((len(dseta), self.N), dtype=np.complex_) # 2D array to store normalized fields
        amu[0, :] = ifft(np.sqrt(2 * self.g / self.kappa) * Amu0) # Store normalized initial field

        # Iterate over all detuning values
        for i in tqdm(range(len(dseta) - 1)):
            dseta_curr = dseta[i] # Current detuning value
            y0 = amu[i, :] # Initial conditions to solve LLE
            if i == 0:
                def LLE(tau, y):
                    amu_ = y
                    damu_ = -(1 + 1j * (dseta_curr + (dseta_step * (tau / tau_step)) + self.dint)) * amu_ + 1j * ifft(np.abs(fft(amu_))**2 * fft(amu_)) + self.f # LLE
                    dy = damu_
                    return dy 
            solver = solve_ivp(LLE, [0, tau_step], y0) # Solve LLE
            amu[i + 1, :] = solver.y[:, -1] # Store field solution
        
        return dseta, amu
