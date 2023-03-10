classdef Ring
    properties
        c = 299792458; % Light speed [m/s]
        hbar = 1.0545718176461565e-34; % Reduced Planck constant [W*s^2]
        N; n0; n2; FSR; lambda0; kappa; eta; Veff; D2; Pin; Tr; freq0; omega0; g; f; d2; mu; dint;
    end
    methods
        function obj = Ring(ring_parameters)
            % Retrieve ring parameters
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
            obj.f(round(obj.N / 2)) = sqrt((8 * obj.g * obj.eta / obj.kappa^2) * (obj.Pin / (obj.hbar * obj.omega0))); % Normalized pump field
            obj.D2 = 2 * pi * obj.D2; % Second order dispersion [rad/s]
            obj.d2 = (2 / obj.kappa) * obj.D2; % Normalized second order dispersion
            obj.mu = -(obj.N - 1) / 2 : (obj.N - 1) / 2; % Mode numbers relative to the pumped mode
            obj.dint = (obj.d2 / 2) * obj.mu.^2; % Normalized integrated dispersion
        end
        function [dseta, amu, theta] = numerical_simulation(obj, parameters, simulation_options, varargin)
            % Retrieve simulation options
            Effects = simulation_options("Effects"); % None or "Thermal" or "Avoided mode crossings"
            Noise = simulation_options("Noise"); % True or False (White noise)

            % Retrieve simulation parameters
            dseta_start = parameters("dseta_start"); % Normalized detuning start
            dseta_end = parameters('dseta_end'); % Normalized detuning end
            dseta_step = parameters('dseta_step'); % Tuning step
            roundtrips_step = parameters('roundtrips_step'); % Roundtrips per tuning step
            if isKey(parameters, 'Amu0'); Amu0 = parameters('Amu0'); else; Amu0 = NaN; end % Initial field
            if Effects == "Thermal"
                if isKey(parameters, 'theta0'); theta0 = parameters('theta0'); else; theta0 = NaN; end % Normalized initial variation of temperature
                tauT = parameters("tauT"); % Thermal relaxation time [s]
                n2T = parameters("n2T"); % Coefficient of thermal nonlinearity [m^2/W]
            elseif Effects == "Avoided mode crossings"
                [theta0, tauT, n2T] = deal(0, 0.1, 0);
                mode_perturbated = parameters('mode_perturbated'); % Position of the modal crossing
            else
                [theta0, tauT, n2T] = deal(0, 0.1, 0);
            end

            % Retrive forward values if available
            try
                dseta_forward = varargin{1};
                amu_forward = varargin{2};
                theta_forward = varargin{3};
            catch
                [dseta_forward, amu_forward, theta_forward] = deal([], [], []);
            end
 
            % Calculate further parameters
            dseta = dseta_start : dseta_step : dseta_end; % Normalized detuning (vector)
            tau_step = (obj.kappa / 2) * obj.Tr * roundtrips_step; % Normalized time per tuning step
            amu = zeros(length(dseta), obj.N); % 2D array to store normalized fields
            if isempty(amu_forward); amu(1, :) = ifft(sqrt(2 * obj.g / obj.kappa) * Amu0); else; amu(1, :) = amu_forward(find(dseta_forward >= dseta_start, 1), :); end % Store normalized initial field
            theta = zeros(length(dseta), 1); % 1D array to store initial normalized variation of temperature
            if isempty(theta_forward); theta(1) = theta0; else; theta(1) = theta_forward(find(dseta_forward >= dseta_start, 1)); end % Store initial normalized variation of temperature
            if Effects == "Thermal"; aux = 1i; else; aux = 0; end % Auxiliary variable to consider or not the variation of temperature
            if Effects == "Avoided mode crossings"; obj.dint(round(obj.N / 2) + mode_perturbated) = 0; end % Perturbate integrated dispersion

            % Iterate over all detuning values
            if isempty(amu_forward); tuning = "FORWARD"; else; tuning = "BACKWARD"; end
            fprintf("%s TUNING: \n", tuning); % Indicate if forward or backward tuning simulation
            progress = 0; % Variable to store simulation progress
            disp(['Progress: ', num2str(progress), '%']); % Display simulation progress
            for i = 1 : (length(dseta) - 1)
                current_progress = round(i / (length(dseta) - 1) * 100); % Current simulation progress
                if progress < current_progress
                    progress = current_progress; % Update simulation progress
                    disp(['Progress: ', num2str(progress), '%']); % Display simulation progress
                end
                dseta_current = dseta(i); % Current detuning value
                y0 = [amu(i, :), theta(i)]; % Initial conditions to solve LLE
                [~, y] = ode45(@CMEs, [0, tau_step], y0); % Solve LLE
                if Noise == "true"; noise = ifft(sqrt(2 * obj.g / obj.kappa) * (randn(1, obj.N) + 1i * randn(1, obj.N))); else; noise = 0; end % White noise
                amu(i + 1, :) = y(end, 1 : end - 1) + noise; % Store field
                theta(i + 1) = y(end, end); % Store variation of temperature
            end

            function dy = CMEs(tau, y) % Function for LLE and thermal equation
                amu_ = y(1 : end - 1); % Retrieve field
                theta_ = y(end); % Retrieve variation of temperature
                damu_ = -(1 + 1i * (dseta_current + dseta_step * (tau / tau_step) + obj.dint.') - aux * theta_) .* amu_ + 1i * ifft(abs(fft(amu_)).^2 .* fft(amu_)) + obj.f.'; % LLE
                dtheta_ = (2 / obj.kappa / tauT) * (n2T / obj.n2 * sum(abs(amu_).^2) - theta_); %Thermal equation
                dy = [damu_; dtheta_];
            end
        end
    end
end