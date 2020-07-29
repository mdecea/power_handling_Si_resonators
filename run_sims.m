% --------------- Non-linearities in a microring resonator --------------

% The script implements and solves a system of differential equations that
% models the dynamics of a microring resonator under high power conditions,
% in which case nonlinearities such as two photon absorption become
% important.

% The details of the model are described in:
% de Cea, Atabaki, Ram: "Power handling of silicon microring modulators", 
% Optics Express 27, n. 17 (2019). 
% DOI: https://doi.org/10.1364/OE.27.024274

% In this script, we parallelize both input power and wavelength. We can
% do this because the algortithm that computes the initial conditions given
% a wavelength and a power is pretty good.

clear all
close all

warning off

%% ********************* 1. Operational conditions ************************

T0 = 300;  % Operating temperature

L = 7; % Length of the PRBS signal (2^L-1)
n = 10; % Number of times the PRBS signal is applied 8
Vbias = -2.5; % Reverse bias voltage (V) -2.5
Vp = 2; % Peak amplitude of the applied driving signal (V)
f = 1e9; % Data rate of the signal (bps)

sim_name = 'trial';  %Name that the .mat files containing sim data will have

% The parameters for the ring shoul be specified in the script
% get_ring_params

% Sweep over input power and wavelength. Again, every [power,wavelength]
% pair will be computed in parallel.

Pin_v = 2.5e-3; % Input power (W)
% lamL_v = linspace(1550.1, 1550.5, 5)*1e-9; % Laser wavelengths (m)
lamL_v = [1550.25]*1e-9;

 
%% ******* 2. Generate PRBS driving signal from specified parameters ******

% Calculate the end_time for each operational point
% from the above parameters
end_t = (2^L-1)*n/f;
samples_per_bit = 20;
nsamples = (2^L-1)*samples_per_bit; % Total number of samples per prbs trace
sample_t = linspace(0, end_t, nsamples*n); % Sampling times

% Generate the driving signal
% Vapp = (Vbias + Vp*idinput([nsamples, 1, n], 'prbs', [0, (2^L-1)/(nsamples)], [-1, 1]));
Vapp = (Vbias + Vp*LUT_PRBS(L, samples_per_bit, n));

%% ***** 3. Perform the simulations for each operational condition ******

% First we need to generate all the combinations of Pin and wavelength
[Pin_m, lamL_m] = meshgrid(Pin_v, lamL_v);

Pin_m = Pin_m(:);
lamL_m = lamL_m(:);

ER = zeros(1, length(Pin_m));
IL = zeros(1, length(Pin_m));
mu_0 = zeros(1, length(Pin_m));
mu_1 = zeros(1, length(Pin_m));
mu_0_sigma = zeros(1, length(Pin_m));
mu_1_sigma = zeros(1, length(Pin_m));

c = 2.997e8;

% Do each power and wavelength in parallel
parfor j = 1:length(Pin_m)
    
    Pin = Pin_m(j);
    lamL = lamL_m(j);
    Wl = 2*pi*c/lamL;
    
    % *******************************
    % Get the initial state as if we have swept the laser wavelength
    lamL_sweep = linspace(lamL-1e-9, lamL, 100);
    Wl_sweep = 2*pi*c./lamL_sweep;
    
    init_guess = [NaN, NaN];
    
    for k = 1:length(Wl_sweep)
        [y0, ~, success] = get_steady_state(Wl_sweep(k), Pin, T0, Vbias, init_guess);
        init_guess = [y0(1), y0(2)];
        fprintf('Init %d out of %d done \n', k, 100);   
    end
    fprintf('Inital state done \n')
       
    % ************************************
    
    % Solve the ODEs for the operational point

    [t, y] = run_single_sim(Wl, Vbias, T0, Pin, y0, sample_t, Vapp, 4*n);

    % Save the relevant results
    [gamma_0, ~, ~, kappa, ... % loss params
    ~, ~, ~, ~, ~, ~, ... % Optical mode volumes and confinements
    ~, ~, ~, ~, ... % Thermal related stuff, instantaneous value at T+deltaT
    ~, ~, ~, ... % Thermal related stuff, equivalent value (see comments in code)
    ~, ~, ~, ~, ~, ... % Carrier related stuff
    ~, ~, ~, ... % Silicon related parameters
    ~,~, ~, ... % Ring related stuff
    ~, ~] ...  % Electro-optic driving related stuff
    = ring_params(Wl, T0, 0, 0, 0);    
    
    Pout = abs(sqrt(Pin) - 1i*conj(kappa).*(y(:,1)+1i*y(:,2))).^2;
    parsave(strcat('data/', sim_name, '_Pout_for_Pin=', num2str(Pin), '_lam=', num2str(lamL*1e9), '.mat'), Pout, t);
    parsave(strcat('data/', sim_name, '_y_for_Pin=', num2str(Pin), '_lam=', num2str(lamL*1e9), '.mat'), y, t);
    
    % Get metrics
    file_name = strcat('data/', sim_name, '_hists_Pin=', num2str(Pin), '_lam=', num2str(lamL*1e9));
    [ER_s, IL_s, mu_0_s, mu_1_s, mu_0_a, mu_1_a] = analyze_traces(Pout, Pin, y(:, 5), Vbias, Vp, file_name);

    ER(j) = ER_s;
    IL(j) = IL_s
    mu_0(j) = mu_0_s;
    mu_0_sigma(j) = mu_0_a;
    mu_1_sigma(j) = mu_1_a;
    mu_1(j) = mu_1_s;

end

% Finally, save the metric data
save(strcat('data/', sim_name, '_metrics.mat'), 'Pin_m', 'lamL_m', 'ER', 'IL', 'mu_0', 'mu_1', 'f', ...
    'mu_0_sigma', 'mu_1_sigma');


%% ****** Helper functions **********

function parsave(fname, data, t)

  save(fname, 'data', 't')
  
end

function [ER, IL, mu_0, mu_1, mu_0_amp, mu_1_amp] = analyze_traces(Pout, Pin, Vpn, Vbias, Vp, filename)
                   
   valV1 = Pout(Vpn > (Vbias + Vp - 0.001));
   valV0 = Pout(Vpn < (Vbias - Vp + 0.001));
   
   % Get relevant values by taking histograms
   
   % '1' value
   figure()
   h = histogram(valV1, 100, 'Normalization', 'probability');
   hold on
   bins = h.BinEdges;
   bins = (bins(2:end)+bins(1:end-1))/2;
   counts = h.Values;
   [~, ind_1] = max(counts);
   mu_1 = bins(ind_1);
   norm1 = fitdist(valV1, 'Normal');
   plot(bins, pdf(norm1,bins)*max(counts)/max(pdf(norm1,bins)), '--', 'LineWidth', 2)
   mu_1_amp = norm1.sigma;
   
   % '0' value
   h = histogram(valV0, 100, 'Normalization', 'probability');
   title(['Pin = ', num2str(Pin*1e3), ' mW'])
   bins = h.BinEdges;
   bins = (bins(2:end)+bins(1:end-1))/2;
   counts = h.Values;
   [~, ind_0] = max(counts);
   mu_0 = bins(ind_0);
   norm0 = fitdist(valV0, 'Normal');
   plot(bins, pdf(norm0,bins)*max(counts)/max(pdf(norm0,bins)), '--', 'LineWidth', 2)
   mu_0_amp = norm0.sigma;
   
   savefig(strcat(filename, '.fig'));
   close(gcf)

   
   if mu_1 < mu_0
       
      mu_int = mu_0;
      mu_0 = mu_1;
      mu_1 = mu_int;
      
      mu_amp_int = mu_0_amp;
      mu_0_amp = mu_1_amp;
      mu_1_amp = mu_amp_int;
      
   end

   ER = mu_1/mu_0;
   IL = Pin/mu_1;

end




