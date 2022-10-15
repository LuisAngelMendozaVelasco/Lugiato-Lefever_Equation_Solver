from scipy.constants import c, hbar, pi
import numpy as np

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
        self.f = np.sqrt((self.Pin / (hbar * self.omega0)) * (8 * self.g * self.eta / self.kappa**2)) # Normalized pump field
        self.f = np.insert(np.zeros(self.N), int((self.N / 2) + 1), self.f) # Normalized pump field vector
        self.D2 = 2 * pi * self.D2 # Second order dispersion [rad/s]
        self.d2 = (2 / self.kappa) * self.D2 # Normalized second order dispersion
        self.mu = np.arange(-(self.N - 1) / 2, ((self.N - 1) / 2) + 1, 1) # Mode numbers relative to the pumped mode
        self.dint = (self.d2 / 2) * self.mu**2 # Normalized integrated dispersion
        