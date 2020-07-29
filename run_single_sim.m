function [t, y] = run_single_sim(Wl, Vbias, T0, Pin, y0, tapp, vapp, break_times)

% This code runs a single simulation for the dynamics of the ring resonator

% Input parameters
% ----------------
% Wl --> Laser frequency (rad/s)
% Vbias --> Bias voltage (V)
% T0 --> Ambient temperature (K)
% Pin --> Input power (W)
% y0 --> Initial conditions for relevant variables: [real(a), im(a), deltaT, Ntpa, Vbias]
% tapp --> sampling times for the applied voltage
% vapp --> Voltage applied at times tapp
% break_times --> Number of time periods in which we break the whole time
    % series. This is for speeding up the simulation.

% Output parameters
% -----------------
% t --> Time vector
% y --> Relevant values at times t. [real(a), im(a), deltaT, Ntpa, Vbias]


% **********************  Computation ****************************

t = [];
y = [];

% As we use interpolation in the nlodesystem funtion, if we 
% interpolate over a very long vector the interpolation takes a lot of time. 
% To avoid this, we break all the simulation in smaller time periods so the 
% interpolation does not need to be over all the vector, thus reducing
% computation time

nsamples = length(vapp);

params = [Wl, Vbias, Pin, T0];

for k = 1:break_times
    tic;
    V_int = vapp(((nsamples/break_times)*(k-1)+1):(nsamples*k/break_times));
    sample_t_int = tapp(((nsamples/break_times)*(k-1)+1):(nsamples*k/break_times));
    
    % Solve the ODE numerically
    [t_part,y_part] = ode113(@(t,y) odesystem(t, y, params, V_int, sample_t_int), ... 
               sample_t_int, y0, odeset('InitialStep', 1e-18, 'RelTol', 3e-12, ...
               'AbsTol', [1e-15, 1e-15, 1e-6, 100, 1e-3])); % 'OutputFcn',@odeplot; 'MaxStep', 0.1e-12, 

    t = [t; t_part];
    y = [y; y_part];
    y0 = y_part(end,:);
    toc;

end


end