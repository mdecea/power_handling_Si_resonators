function y = equilibrium(x, Wl, Pin, T0)
    
    % Obtains the value of the electric field inside the ring in steady state. 
    % When equilibrium(x) = 0 we are at equilibrium.

    % x(1) --> Real part of the electric field inside the ring
    % x(2) --> Imaginary part of the electric field inside the ring
    
    c = 2.997e8;
    hbar = 1.05e-34;
    
    % Solve for deltaT iteratively until we converge
    deltaT = 0;
    prev_sol = Inf;
    
    num_counts = 0;
    
    while abs(deltaT-prev_sol) > 1 && num_counts < 100
        
        prev_sol = deltaT;
        
        % Get deltaT from the provided x, assuming we are close to T0
        
        [gamma_0, gamma_lin, gamma_rad, kappa, ... % loss params
        Vring, V_TPA, V_FCA, gamma_disk, gamma_FCA, gamma_TPA, ... % Optical mode volumes and confinements
        ~, rho_Si, ~, ~, ... % Thermal related stuff, instantaneous value at T+deltaT
        gamma_th_eq, Cp_Si_eq, dnSidT_eq, ... % Thermal related stuff, equivalent value (see comments in code)
        gamma_r, alfa_p, alfa_n, dnSidNn, dnSidNp, ... % Carrier related stuff
        beta_Si, n_Si, ng, ... % Silicon related parameters
        Wl, W0,  deltaW0_cold, ... % Ring related stuff
        dW0dV, tau] ...  % Electro-optic driving related stuff
        = ring_params(Wl, T0, deltaT, 0, 0);
        
        % For equilibrium considerations, we care abot the thermal
        % equivalent values gamma_th_eq, Cp_Si_eq, dnSidT_eq
        H = (gamma_FCA*beta_Si*c^2)/(2*hbar*(Wl)*ng^2*V_FCA^2*gamma_r);
        G = (gamma_disk/(rho_Si*Cp_Si_eq*Vring*gamma_th_eq));
        E = (gamma_TPA*beta_Si*c^2/(V_TPA*ng^2));
        F = (alfa_p+alfa_n)*c*H/ng;
        
        a_re = x(1)*1e-16;
        a_im = x(2)*1e-16;
        u = a_re.^2 + a_im.^2;
        deltaT = G*(gamma_lin + E*u + F*u^2)*u;
        
        %abs(deltaT-prev_sol)
        %pause
        
        num_counts = num_counts+1;
        
    end

    lam = (gamma_rad + gamma_0 + gamma_lin + E*u + F*u^2)/2;

    B = deltaW0_cold - (-W0/n_Si)*(dnSidT_eq*G*(gamma_lin+E*u+F*u^2)*u + dnSidNn*(H*u^2*1e-6)^1.05 + dnSidNp*(H*u^2*1e-6)^0.8);

    y(1) = - lam*a_re - B*a_im;
    y(2) = -lam*a_im + B*a_re - kappa*sqrt(Pin);

end