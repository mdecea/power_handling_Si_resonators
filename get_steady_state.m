function [y0, deltaW0, success] = get_steady_state(Wl, Pin, T0, Vbias, init_guess)
    
    % Gets the steady state solution for the system of nonlinear equations
    
    % Inputs
    % ------
    % % Wl --> Laser frequency (rad/s)
    % Pin --> Input power (W)
    % T0 --> Ambient temperature (K)
    % Vbias --> Bias voltage (V)
    % init_guess (optional) --> If specified, uses init_guess as the 
    %       initial guess for the solution of the nonlinear system
   
    c = 2.997e8;
    hbar = 1.05e-34;
    
    if nargin == 5 && ~ logical(sum((isnan(init_guess))))
        a_steady = init_guess;
    else
        [gamma_0, gamma_lin, gamma_rad, kappa, ... % loss params
        Vring, V_TPA, V_FCA, gamma_disk, gamma_FCA, gamma_TPA, ... % Optical mode volumes and confinements
        gamma_th, rho_Si, Cp_Si, dnSidT, ... % Thermal related stuff, instantaneous value at T+deltaT
        gamma_th_eq, Cp_Si_eq, dnSidT_eq, ... % Thermal related stuff, equivalent value (see comments in code)
        gamma_r, alfa_p, alfa_n, dnSidNn, dnSidNp, ... % Carrier related stuff
        beta_Si, n_Si, ng, ... % Silicon related parameters
        Wl, W0,  deltaW0_cold, ... % Ring related stuff
        dW0dV, tau] ...  % Electro-optic driving related stuff
        = ring_params(Wl, T0, 0, 0, 0);
           
        a_steady = 1i*kappa*sqrt(Pin)/((gamma_rad + gamma_0 + gamma_lin)/2+1i*deltaW0_cold); % Initial guess for the solver
        % a_steady = [0, 0];
        
        if deltaW0_cold == 0
           Wl = Wl + 1e12;
        end
    end
    

    % Solve the trascendental equation
    opt = optimoptions('fsolve', 'TolFun', 1e-10, 'StepTolerance', 1e-12, 'Display', 'off');
    [a_steady, ~, exitflag, ~] = fsolve(@(x) equilibrium(x, Wl, Pin, T0), [real(a_steady)*1e16 imag(a_steady)*1e16], opt);
    
    a_steady = a_steady*1e-16;
    
    % Refine the solution iteratively to be self consistent with
    % deltaT.
    deltaT = 0;
    prev_sol = Inf;
    num_counts = 0;
    
    while abs(deltaT-prev_sol) > 0.5 && num_counts < 100
        
        prev_sol = deltaT;
        
        [gamma_0, gamma_lin, gamma_rad, kappa, ... % loss params
        Vring, V_TPA, V_FCA, gamma_disk, gamma_FCA, gamma_TPA, ... % Optical mode volumes and confinements
        gamma_th, rho_Si, Cp_Si, dnSidT, ... % Thermal related stuff, instantaneous value at T+deltaT
        gamma_th_eq, Cp_Si_eq, dnSidT_eq, ... % Thermal related stuff, equivalent value (see comments in code)
        gamma_r, alfa_p, alfa_n, dnSidNn, dnSidNp, ... % Carrier related stuff
        beta_Si, n_Si, ng, ... % Silicon related parameters
        Wl, W0,  deltaW0_cold, ... % Ring related stuff
        dW0dV, tau] ...  % Electro-optic driving related stuff
        = ring_params(Wl, T0, deltaT, 0, 0);
    
    
        % This is an equilibrium calcualtion, so we care about thermal
        % equivalent values
        H = (gamma_FCA*beta_Si*c^2)/(2*hbar*(Wl)*ng^2*V_FCA^2*gamma_r);
        G = (gamma_disk/(rho_Si*Cp_Si_eq*Vring*gamma_th_eq));
        E = (gamma_TPA*beta_Si*c^2/(V_TPA*ng^2));
        F = (alfa_p+alfa_n)*c*H/ng;
        
        a_re = a_steady(1);
        a_im = a_steady(2);
        u = a_re.^2 + a_im.^2;
        deltaT = G*(gamma_lin + E*u + F*u^2)*u;
        
        %abs(deltaT-prev_sol)
        num_counts = num_counts+1;
        
    end
     
    Ntpa_steady = H*u^2;
    deltaT_steady = G*(gamma_lin + E*u + F*u^2)*u;
    
    deltaW0 = deltaW0_cold - (-W0/n_Si)*(dnSidT_eq*G*(gamma_lin+E*u+F*u^2)*u ...
        + dnSidNn*(H*u^2*1e-6)^1.05 + dnSidNp*(H*u^2*1e-6)^0.8);
        % Detuning between laser wavelenght and resonance wavelength = 
        % = Wl - W0
    
    y0 = [a_steady(1), a_steady(2), deltaT_steady, Ntpa_steady, Vbias];
    
    if exitflag < 1
        success = 0;
    else
        success = 1;
    end
    
    %s_out = sqrt(Pin) - 1i*kappa*(a_steady(1)+1i*a_steady(2));
    %P_out = abs(s_out)^2

end