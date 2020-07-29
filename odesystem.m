function dydt = odesystem(t, y, params, V_app, t_app)

% This function returns the derivative with respect to time of relevant
% variables

%% *********** Get parameters ************
Wl = params(1);
Vbias = params(2);
Pin = params(3);
T0 = params(4);

[gamma_0, gamma_lin, gamma_rad, kappa, ... % loss params
    Vring, V_TPA, V_FCA, gamma_disk, gamma_FCA, gamma_TPA, ... % Optical mode volumes and confinements
    gamma_th, rho_Si, Cp_Si, dnSidT, ... % Thermal related stuff, instantaneous value at T+deltaT
    gamma_th_eq, Cp_Si_eq, dnSidT_eq, ... % Thermal related stuff, equivalent value (see comments in code)
    gamma_r, alfa_p, alfa_n, dnSidNn, dnSidNp, ... % Carrier related stuff
    beta_Si, n_Si, ng, ... % Silicon related parameters
    Wl, W0,  deltaW0_cold, ... % Ring related stuff
    dW0dV, tau] ...  % Electro-optic driving related stuff
    = ring_params(Wl, T0, y(3), y(4), y(5));

% % Voltage applied 
% Obtain the voltage applied to the electrodes at time t by nearest
% neighbor interpolation (which is good enough since we are applying a
% square wave). To do it, we implement an ordered insertion algorithm.

lo = 1;
hi = length(t_app);
while (lo < hi)
    mid = floor((lo+hi)/2);
    if t < t_app(mid)
        hi = mid;
    else
        lo = mid + 1;
    end
end

V_t = V_app(lo); % Voltage applied at time t

%% *********** Constants  ***********

hbar = 1.05e-34;  % Reduced Planck constant (J*s)
c = 2.997e8; % Speed of light (m/s)


%% ***********  Relevant variables *********

% For this, we use the equivalent thermal constants

deltaW0_mod = dW0dV*(y(5)-Vbias); % Dispersion due to applied voltage

deltaW0 = (((-1/n_Si)*(dnSidT_eq*y(3) + ...
    dnSidNn*(y(4)*1e-6)^1.05 + ...
    dnSidNp*(y(4)*1e-6)^0.8 )) + ...
    deltaW0_mod/W0)*W0;  % Overall dispersion (W0(t) - W0_cold)

U = abs(y(1)+1i*y(2)).^2;  % Energy in the ring

TPA_loss = (gamma_TPA*beta_Si*c^2/(V_TPA*ng^2))*U;  % TPA optical loss
FCA_loss = (alfa_p*y(4) + alfa_n*y(4))*c/ng; % FCA optical loss
gamma = gamma_rad + gamma_0 + gamma_lin + TPA_loss + FCA_loss;  % Total loss rate

Pabs = (gamma_lin + TPA_loss + FCA_loss)*U;  % Absorbed power

G = (gamma_FCA*beta_Si*c^2/(2*hbar*(Wl)*ng^2*V_FCA^2))*U.^2;  % TPA generation


%% ********** Derivatives as given by the model **********

% Calculate the derivatives as given by the model
da = (-gamma/2 + 1i*(deltaW0_cold - deltaW0)).*(y(1)+1i*y(2)) - 1i*kappa*(sqrt(Pin));
dT = -gamma_th*y(3) + (gamma_disk/(rho_Si*Cp_Si*Vring))*Pabs;
dNtpa = - gamma_r*y(4) + G;
dV = -y(5)/tau + V_t/tau;  % Voltage at the junction


%% ******* Construct vector ******

dydt = [real(da); imag(da); dT; dNtpa; dV];


end



