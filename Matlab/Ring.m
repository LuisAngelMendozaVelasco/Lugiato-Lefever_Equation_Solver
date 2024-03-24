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
            wb = waitbar(0, '0.0%', 'Name', 'Simulation progress', 'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
            setappdata(wb, 'canceling', 0);
            fprintf("In progress...\n");
            steps = length(dseta) - 1;

            for i = 1 : steps
                % Check for clicked Cancel button
                if getappdata(wb, 'canceling')
                    fprintf("Canceled!\n");
                    delete(wb);
                    break
                end

                current_progress = round(i / steps * 100, 1); % Current simulation progress

                % Update waitbar and message
                if progress < current_progress
                    progress = current_progress; % Update simulation progress
                    waitbar(i/steps, wb, strcat(sprintf('%.1f', progress), "%"))
                end

                dseta_current = dseta(i); % Current detuning value
                y0 = [amu(i, :), theta(i)]; % Initial conditions to solve LLE
                [~, y] = ode78(@LLE, [0, tau_step], y0); % Solve LLE
                if Noise == "true"; noise = ifft(sqrt(2 * obj.g / obj.kappa) * (randn(1, obj.N) + 1i * randn(1, obj.N))); else; noise = 0; end % White noise
                amu(i + 1, :) = y(end, 1 : end - 1) + noise; % Store field
                theta(i + 1) = y(end, end); % Store variation of temperature

                if i == steps
                    fprintf("Completed!\n");
                    delete(wb);
                end
            end

            function dy = LLE(tau, y) % Function for LLE and thermal equation
                amu_ = y(1 : end - 1); % Retrieve field
                theta_ = y(end); % Retrieve variation of temperature
                damu_ = -(1 + 1i * (dseta_current + dseta_step * (tau / tau_step) + obj.dint.') - aux * theta_) .* amu_ + 1i * ifft(abs(fft(amu_)).^2 .* fft(amu_)) + obj.f.'; % LLE
                dtheta_ = (2 / obj.kappa / tauT) * (n2T / obj.n2 * sum(abs(amu_).^2) - theta_); %Thermal equation
                dy = [damu_; dtheta_];
            end
        end


        function plot_results(obj, dseta_forward, amu_forward, dseta_snap, varargin)
             % Retrive backward values if available
            try
                dseta_backward = varargin{1};
                amu_backward = varargin{2};
            catch
                [dseta_backward, amu_backward] = deal([], []);
            end
 
            figure();
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); %Set figure size

            % Plot average intracavity intensity
            avg_intr_int_fwrd = sum(abs(amu_forward).^2, 2); % Average intracavity intensity of forward tuning [a. u.]
            if ~isempty(amu_backward); avg_intr_int_bwrd = sum(abs(amu_backward).^2, 2); else; avg_intr_int_bwrd = []; end % Average intracavity intensity of backward tuning [a. u.]
            
            subplot(3, 2, 1:2);
            plot(dseta_forward, avg_intr_int_fwrd, Color="b", LineWidth=2);
            hold on;

            if ~isempty(avg_intr_int_bwrd)
                plot(dseta_backward, avg_intr_int_bwrd, Color="r", LineWidth=2);
                xline(dseta_snap, Color="k", LineStyle="--", LineWidth=2);
                legend({'Forward tuning', 'Backward tuning', ['\zeta_0 = ', num2str(dseta_snap)]}, Location="northwest");
            else
                xline(dseta_snap, Color="k", LineStyle="--", LineWidth=2);
                legend({'Forward tuning', ['\zeta_0 = ', num2str(dseta_snap)]}, Location="northwest");
            end

            xlabel('Normalized detuning, $\zeta$', 'interpreter', 'latex');
            ylabel({'Average intracavity'; 'intensity [a. u.]'}); 
            xlim([min(dseta_forward), max(dseta_forward)]);
            ylim([0, max(avg_intr_int_fwrd) * 1.05]);
            grid minor;
            set(gca, 'fontsize', 16);
            
            % Plot forward optical spectrum
            if ~isempty(amu_backward); subplot_ = subplot(3, 2, 3); else; subplot_ = subplot(3, 2, 3:4); end
            plot_optical_spectrum(dseta_forward, amu_forward, dseta_snap, subplot_, "b");

            % Plot forward intracavity power
            if ~isempty(amu_backward); subplot_ = subplot(3, 2, 5); else; subplot_ = subplot(3, 2, 5:6); end
            plot_intracavity_power(dseta_forward, amu_forward, dseta_snap, subplot_, "b");

            if ~isempty(amu_backward)
                % Plot backward optical spectrum
                subplot_ = subplot(3, 2, 4);
                plot_optical_spectrum(dseta_backward, amu_backward, dseta_snap, subplot_, "r");

                % Plot backward intracavity power
                subplot_ = subplot(3, 2, 6);
                plot_intracavity_power(dseta_backward, amu_backward, dseta_snap, subplot_, "r");
            end

            function plot_optical_spectrum(dseta, amu, dseta_snap, subplot_, line_color)
                [~, index] = min(abs(dseta - dseta_snap)); % Index of dseta_snap value in dseta
                omegamu = obj.omega0 + 2 * pi * obj.FSR * obj.mu + (2 * pi * obj.D2 * obj.mu.^2 / 2); % Resonance frequencies [rad/s]
                freqmu = omegamu / (2 * pi); % Resonance frequencies [Hz]
                fin = max(obj.f) * ones(1, obj.N); % Normalized pump field (vector)
                fout = obj.f - 2 * obj.eta * abs(amu(index, :)); % Normalized output field
                optical_spectrum = (abs(fout).^2 ./ abs(fin).^2) * obj.Pin; % Optical spectrum [W]
                optical_spectrum = 10 * log10(1000 * optical_spectrum); % Optical spectrum [dBm]
                
                subplot_;
                stem(freqmu * 1e-12, optical_spectrum, 'Marker', 'none', BaseValue=-30, Color=line_color)
                xlabel('Resonance frequencies [THz], $f_\mu$', 'interpreter', 'latex');
                ylabel({'Optical'; 'spectrum [dBm]'});
                xlim([freqmu(1) * 1e-12, freqmu(end) * 1e-12]);
                ylim([-30, 10 * log10(1000 * obj.Pin) * 1.05]);
                grid minor;
                set(gca, 'fontsize', 16);
            end

            function plot_intracavity_power(dseta, amu, dseta_snap, subplot_, line_color)
                [~, index] = min(abs(dseta - dseta_snap)); % Index of dseta_snap value in dseta
                ring_circumference = linspace(-pi, pi, 4096); % Ring circumference from -pi to pi
                temp = zeros(1, 4096);
                temp(obj.mu + 4096 / 2 + 1) = amu(index, :); % Improve resolution considering more sampling points
                intracavity_power = abs(fft(temp)).^2; % Intracavity power [a. u.]

                subplot_;
                plot(ring_circumference, intracavity_power, Color=line_color, LineWidth=2);
                xlabel('Ring circumference [rad], $\phi$', 'interpreter', 'latex');
                ylabel({'Intracavity'; 'power [a. u]'});
                xlim([-pi, pi]);
                ylim([0, max(intracavity_power) * 1.05])
                xticks([-pi, -pi/2, 0, pi/2, pi]);
                xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
                grid minor;
                set(gca, 'fontsize', 16);
            end
        end
    end
end