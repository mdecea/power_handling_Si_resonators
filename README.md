# Power handling of Si microring resonators

This repository contains the matlab code to solve the time-domain nonlinear model describing nonlinear effects arising in silicon microring modulators.

A detailed description of the model and analysis results are described in the following publication:
Marc de Cea, Amir H. Atabaki, and Rajeev J. Ram, "Power handling of silicon microring modulators," Opt. Express 27, 24274-24285 (2019) Marc de Cea, Amir H. Atabaki, and Rajeev J. Ram, "Power handling of silicon microring modulators," Opt. Express 27, 24274-24285 (2019)

Each function and script is thoroughly commented, so it should be easy to understand the different pieces of the implementation.

The user is suggested to start by looking at the script "run_sims.m". This is the master script. An essential function that the user needs to implement to adapt the code to the specific device is "ring_params.m": this function returns the different parameters necessary for the model to be solved (see Table 1 in the refrenced paper).
