classdef Ring
    properties
        c = 299792458; % Light speed [m/s]
        hbar = 1.0545718176461565e-34; % Reduced Planck constant [W*s^2]
        N; n0; n2; FSR; lambda0; kappa; eta; Veff; D2; Pin; Tr; freq0; omega0; g; f; d2; mu; dint;
    end
    methods
        function obj = Ring(ring_parameters)
            obj.N = ring_parameters("N"); % Number of modes. It must be odd!
            obj.n0 = ring_parameters("n0"); % Refractive index
            obj.n2 = ring_parameters("n2"); % Nonlinear reftactive index [m^2/W]
            obj.FSR = ring_parameters("FSR"); % Free Spectral Range [Hz]
            obj.lambda0 = ring_parameters("lambda0"); % CW pump wavelength [m]
            obj.kappa = ring_parameters("kappa"); % Optical linewidth [Hz]
            obj.eta = ring_parameters("eta"); % Coupling efficiency
            obj.Veff = ring_parameters("Veff"); % Effective mode volume [m^3]
            obj.D2 = ring_parameters("D2"); % Second order dispersion [Hz]
            obj.Pin = ring_parameters("Pin"); % Pump power [W]

            % Calculate further parameters
            obj.Tr = 1 / obj.FSR; % Roundtrip time [s]
            obj.freq0 = obj.c / obj.lambda0; % CW pump frequency [Hz]
            obj.omega0 = 2 * pi * obj.freq0; % CW pump angular frequency [rad/s]
            obj.g = obj.hbar * obj.omega0^2 * obj.c * obj.n2 / (obj.n0^2 * obj.Veff); % Nonlinear coupling coefficient
            obj.kappa = 2 * pi * obj.kappa; % Optical linewidth [rad/s]
            obj.f = zeros(1, obj.N); % Normalized pump field vector
            obj.f(floor(obj.N / 2) + 1) = sqrt((8 * obj.g * obj.eta / obj.kappa^2) * (obj.Pin / (obj.hbar * obj.omega0))); % Normalized pump field
            obj.D2 = 2 * pi * obj.D2; % Second order dispersion [rad/s]
            obj.d2 = (2 / obj.kappa) * obj.D2; % Normalized second order dispersion
            obj.mu = -(obj.N - 1) / 2 : (obj.N - 1) / 2; % Mode numbers relative to the pumped mode
            obj.dint = (obj.d2 / 2) * obj.mu.^2; % Normalized integrated dispersion
        end
    end
end