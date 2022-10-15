from scipy.constants import c, hbar

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