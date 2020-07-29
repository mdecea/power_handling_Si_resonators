function [gamma_0, gamma_lin, gamma_rad, kappa, ... % loss params
    Vring, V_TPA, V_FCA, gamma_disk, gamma_FCA, gamma_TPA, ... % Optical mode volumes and confinements
    gamma_th, rho_Si, Cp_Si, dnSidT, ... % Thermal related stuff, instantaneous value at T+deltaT
    gamma_th_eq, Cp_Si_eq, dnSidT_eq, ... % Thermal related stuff, equivalent value (see comments in code)
    gamma_r, alfa_p, alfa_n, dnSidNn, dnSidNp, ... % Carrier related stuff
    beta_Si, n_Si, ng, ... % Silicon related parameters
    Wl, W0,  deltaW0_cold, ... % Ring related stuff
    dW0dV, tau] ...  % Electro-optic driving related stuff
    = ring_params(Wl, T0, deltaT, deltaN, Vpn)

% This function returns the relevant parameters for the nonlinear model.

% These are all listed in Table 1 of the paper 
% (de Cea, Atabaki, Ram: "Power handling of silicon microring modulators", 
% Optics Express 27, n. 17 (2019) )

% Some parameters are fixed, others depend on temperature and/or free
% carrier concentration.

% Inputs
% ------
% Wl --> Laser freq
% T0 --> Ambient Temperature
% deltaT --> Temperature difference from ambient T (K)
% deltaN --> Free carrier concentration with respect to equilibrium (m^-3)
% Vpn --> Voltage in the pn junction (V)

% Outputs
% -------
% All parameters

%% ********* 0. Constants ***************************

q = 1.6021766208e-19; % Elementary charge (Coulomb)
k = 1.38064852e-23; % Boltzmann constant (J/K)
c = 2.997e8;
e0 = 8.854187817e-12; % Vacuum permitivity (F/m)
hbar = 1.05e-34;

%% ********* 1. Ring parameters ***************************

n_Si = 3.485; % Si index of refraction
ng = 2.769; % Group index

% Ring parameters
lam0 = 1550e-9;
W0 = 2*pi*c/lam0;
deltaW0_cold = Wl-W0;

% Optical mode parameters
gamma_TPA = 0.90967; % Confinement factors (calculated with FD3D mode solver)
gamma_FCA = 0.9621118; % Confinement factors (calculated with FD3D mode solver)
gamma_disk = 0.6515; % Confinement factors (calculated with FD3D mode solver)
V_FCA = 8.67e-18; % m^3
V_TPA = 10.62978e-18; % m^3
Vring = 0; % Ring volume (m^3) (not specified becasue of Non Disclosure)

if Vring == 0
    fprintf('You need to specify the ring volume. Not included because of NDA')
end

rho_Si = 2.33e6; % Density of Si (g/m^3).
                % Source: "Calculation of density and heat capacity of 
                % silicon by molecular dynamics simulation"  
                % R Kojima Endo, Y Fujihara, M Susa
 
alfa_n = 8.5e-22; % FCA coefficient for electrons (Soref) (m^2)
alfa_p = 6e-22; % FCA coefficients for holes (Soref) (m^2)

   
% Parameters taken from the resonance
gamma_0 = 3e10; % losses of the fundamental mode due to coupling to waveguide
gamma_lin = 1.9e10; % Linear absorption coefficient
gamma_rad = 1.1e10; % scattering losses
kappa = sqrt(gamma_0); % Coupling from waveguide to resonator (sqrt(Hz))
gamma_th = 2.4e6; % Thermal decay rate
gamma_r = 10e9; % Recombination decay rate. If fw bias, this has to change!!

% Literature values
dnSidT = 1.86e-4; % thermal derivative of the Si index (K^-1)
dnSidNp = -8.5e-18; % plasma dispersion effect for holes (cm^3)
dnSidNn = -8.8e-22; % plasma dispersion effect for electrons (cm^3)
beta_Si = 8.4e-12; % TPA sthrength (m/W)
Cp_Si = 0.7145; % Specific heat of Silicon (J/(g*K))

dlamdV = -20e-12;  % modulation efficiency (20 pm/V)
dW0dV = 2*pi*c/(lam0+dlamdV) - 2*pi*c/lam0; % Change in the resonnance wavelength due to modulation

% Since at RT the values are not T dependent, the equivalent thermal
% values are the same as the instantanoeus thermal values.
% If these values were dependent on T, the variables xx_eq should be 
% calculated as follows: 
% xx_eq = [integral(xx(T) dT) between T0 and deltaT]/deltaT
gamma_th_eq =  gamma_th;
Cp_Si_eq = Cp_Si;
dnSidT_eq = dnSidT;
    
tau = 1/(2*pi*gamma_r); % RC time constant of the system

end
