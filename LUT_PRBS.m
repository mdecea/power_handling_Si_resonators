function pattern = LUT_PRBS(N, samples_per_bit, num_repeats)

% Generates a PRBS signal of length 2^N-1. Since we have not succeeeded in
% getting the same results as matlab, we will just cheat and copy the
% sequences. 

% Inputs
% ------
% N --> Order of the PRBS (L = 2^N-1)
% samples_per_bit --> Number of samples per each bit
% num_repeats --> Number of times the PRBS pattern is repeated

% Outputs
% -------
% pattern --> The pattern, of total length
% (2^N-1)*samples_per_bit*num_repeats. The signa toggles between 1 and -1.

% ************ Select the correct seed *************
switch N
    
    case 2
    output_pattern = [0 0 1];
    
    case 3
    output_pattern = [0 0 0 1 0 1 1];
    
    case 4
    output_pattern = [0 0 0 0 1 0 1 0 0 1 1 0 1 1 1];
    
    case 5
    output_pattern = [0 0 0 0 0 1 1 0 0 1 0 1 1 0 1 1 1 1 0 1 0 1 0 0 0 1 0 0 1 1 1];
    
    case 6
    output_pattern = [0 0 0 0 0 0 1 0 1 0 1 0 0 1 1 0 0 1 0 0 0 1 0 0 1 0 1 1 0 1 1 0 0 0 1 1 1 0 1 0 0 0 0 1 1 0 1 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1];

    case 7
    output_pattern = [0 0 0 0 0 0 0 1 0 1 0 1 0 1 1 0 0 1 1 0 0 0 1 0 0 ...
        0 1 0 1 1 0 1 0 0 1 1 1 0 0 1 0 0 0 0 1 0 0 1 0 1 0 0 1 0 0 1 1 ...
        0 1 1 0 1 1 1 0 0 0 1 1 1 1 0 1 0 0 0 0 0 1 1 0 1 0 1 0 0 0 1 1 ...
        0 0 1 0 1 1 1 0 1 1 0 0 0 0 1 1 1 0 1 0 1 1 1 1 0 0 1 1 1 1 1 0 ...
        1 1 1 1 1 1];
    
    case 8
    output_pattern = [0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 1 0 0 0 0 1 1 1 0 0 1 0 1 0 0 0 1 1 0 1 1 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0 1 0 1 1 1 1 1 0 1 1 1 1 1 1 0 1 0 0 1 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1 1 0 1 0 0 0 1 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1 1 0 1 1 0 1 1 1 0 1 0 1 1 1 1 0 1 1 0 0 1 0 1 1 1 0 0 0 0 0 1 0 1 0 0 1 0 1 1 0 1 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 1 0 0 1 0 1 0 1 1 1 0 1 1 0 1 0 0 0 0 1 1 0 1 0 0 1 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 1 1 1 1 1 0 0 0 1 0 0 1 1 1 1 0 1 0 1 0 0 1 1 0 1 1 0 0 0 1 1 0 0 0 1 0 1 0 1];
    
    case 9
    output_pattern = [0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 1 1 1 1 0 1 0 0 1 1 0 0 1 0 0 1 0 0 0 0 1 0 1 1 1 1 0 0 0 1 1 0 0 1 1 1 1 0 1 1 0 1 1 1 0 1 0 1 0 0 0 1 0 1 0 0 0 0 1 1 0 1 1 0 1 0 0 0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 1 0 0 0 1 0 1 1 0 0 0 0 1 0 1 0 1 1 0 1 0 1 1 1 1 1 1 0 1 0 1 0 1 0 1 0 0 0 0 0 1 0 1 0 0 1 0 1 1 1 1 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 1 1 1 1 0 1 0 0 0 1 0 0 0 0 0 1 1 1 0 0 0 0 1 1 0 0 1 0 1 1 0 0 1 0 1 0 0 0 1 1 1 0 0 1 0 1 1 1 0 1 0 0 0 0 0 0 0 1 0 1 1 0 1 0 0 1 1 1 0 1 0 1 1 0 0 1 1 1 0 0 1 1 1 1 1 1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 1 0 1 1 0 1 1 0 1 1 0 0 1 0 0 0 0 0 0 1 1 0 1 0 0 1 0 1 0 1 1 1 1 0 1 0 1 1 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0 0 0 0 1 0 0 1 1 1 1 0 0 1 0 1 0 1 0 1 1 0 0 0 1 1 0 1 1 1 1 0 0 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 1 0 1 1 1 0 0 0 1 0 1 0 1 0 0 1 0 0 1 1 1 0 0 0 1 1 1 0 1 1 0 1 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 1 1 1 0 1 0 0 1 0 0 0 1 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 0 1 1 0 0 0 1 0 1 1 1 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1];
    
    case 10
    output_pattern = [0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0 1 0 0 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 0 0 1 0 1 1 0 1 0 1 0 1 1 1 1 1 0 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 1 0 0 0 0 1 0 1 1 1 1 0 0 0 1 0 1 1 0 1 1 1 0 0 1 1 0 1 0 0 1 0 1 0 0 1 1 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 1 0 0 0 0 0 0 1 1 0 1 0 1 0 1 0 1 1 0 0 1 1 0 0 1 1 0 1 0 1 1 0 0 0 0 0 1 0 1 1 0 0 0 1 1 1 1 0 1 1 1 0 0 1 0 0 1 1 0 1 1 1 0 1 0 1 1 0 0 1 0 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 0 1 1 0 0 0 1 0 0 0 1 0 0 0 1 1 0 0 0 1 0 1 0 1 1 0 0 0 1 0 1 1 1 1 1 0 0 0 0 1 0 0 1 0 0 0 1 1 1 1 0 0 1 1 1 0 1 1 0 1 0 1 1 0 1 0 0 1 1 0 0 1 0 1 1 1 0 1 1 1 0 1 0 0 1 0 1 1 0 1 0 0 0 1 0 1 1 0 0 1 1 1 0 1 0 0 1 1 1 1 1 1 0 1 0 1 1 0 1 1 0 1 0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 1 1 1 0 0 1 0 0 0 1 0 0 1 1 1 1 0 0 0 0 1 1 0 1 1 0 0 0 1 1 0 1 0 0 1 1 1 0 1 1 1 1 0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 0 1 0 1 1 0 1 0 1 1 1 1 0 1 1 1 1 0 1 1 0 1 0 0 1 0 0 0 0 0 1 0 1 0 0 0 1 1 1 0 1 0 0 0 1 1 0 1 1 1 1 0 0 0 0 0 1 0 0 1 0 1 0 1 0 1 1 1 0 1 0 0 0 0 1 0 0 1 1 0 0 0 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 0 1 0 1 0 0 1 1 0 1 0 0 0 0 1 1 0 1 0 0 0 1 1 1 1 1 0 1 0 1 0 0 1 0 0 1 1 0 0 1 1 1 1 0 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 0 1 0 0 0 0 0 0 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 0 1 1 0 1 0 1 1 1 0 0 1 0 1 1 1 1 1 1 0 0 1 1 0 1 1 0 1 1 1 0 1 1 1 1 1 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 1 1 0 0 1 0 1 0 1 0 0 1 1 1 1 0 1 0 0 0 1 0 0 1 0 1 1 1 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 1 0 0 0 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 0 1 0 0 1 0 0 1 1 1 0 1 0 1 1 1 0 1 1 0 0 1 1 0 1 1 1 1 1 0 0 1 0 1 1 0 1 1 0 0 0 0 1 0 0 0 0 0 1 1 1 0 1 0 1 0 1 0 0 1 0 1 1 1 1 0 1 0 1 1 1 1 1 1 1 0 1 0 0 1 0 0 1 0 0 0 0 1 1 0 0 0 0 1 1 1 0 1 1 1 0 0 0 0 0 0 1 0 0 1 1 1 0 0 0 1 0 1 0 0 1 0 1 0 1 1 1 1 0 0 1 1 0 0 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0 0 1 0 0 1 0 0 1 0 1 0 0 0 1 0 1 0 0 0 0 1 1 1 1 0 1 0 1 0 1 1 0 1 1 1 1 0 1 0 0 1 1 0 1 1 0 0 1 1 1 1 1 0 1 1 1 0 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 0 1 1 0 0 1 0 1 1 0 0 1 0 1 0 0 0 0 0 1 1 0 0 1 1 1 0 0 0 0 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 0 0 1 1 1 1 0 0 0 1 1 1 1 1 1 1];
    
    case 11
    output_pattern = [0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 1 0 0 1 1 1 1 1 1 0 1 0 0 0 1 1 1 1 0 1 1 0 1 0 0 1 1 0 1 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 1 1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 0 0 1 1 1 0 1 0 1 1 1 0 0 1 0 1 1 1 0 1 0 0 0 1 1 0 1 0 1 1 0 1 0 0 0 1 1 1 0 0 1 1 0 1 0 0 1 0 0 0 0 0 1 1 0 0 1 0 1 1 1 0 0 0 0 0 1 1 0 1 0 0 1 1 1 0 0 0 1 1 0 0 0 1 0 0 1 0 0 0 0 1 0 1 0 0 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 0 1 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 1 0 0 1 0 1 0 1 1 1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 0 0 1 0 0 0 0 1 0 1 1 0 1 0 1 1 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 1 1 0 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1 1 1 0 0 0 0 1 0 1 0 1 0 0 1 1 0 1 1 1 1 1 0 0 0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 1 0 0 0 0 0 1 0 0 1 1 0 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 0 1 1 0 1 1 1 1 0 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1 1 0 1 1 1 0 1 1 0 1 0 1 0 1 0 1 0 0 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 1 0 0 0 1 1 0 0 1 0 1 0 1 0 0 0 0 0 1 1 1 1 1 0 1 1 1 0 0 1 1 1 0 1 0 1 0 0 0 0 1 0 1 1 1 1 0 1 1 0 1 1 0 1 1 0 1 0 0 1 0 0 1 0 0 1 1 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 1 0 1 1 1 1 0 1 0 0 1 1 0 1 1 0 1 1 0 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 1 1 0 0 1 1 1 0 1 0 0 0 0 0 0 1 0 1 1 0 1 1 1 1 0 1 1 0 0 1 0 1 1 0 1 0 0 0 0 1 1 0 0 1 1 0 1 1 0 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 0 1 0 1 0 0 1 1 1 0 1 1 1 1 0 0 0 1 0 1 0 1 1 0 0 1 0 1 1 1 1 0 0 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 1 0 1 0 1 1 1 1 1 0 1 1 1 1 0 1 1 1 0 1 0 1 1 0 1 0 1 0 1 1 1 0 0 1 1 1 1 1 0 1 0 0 0 0 1 1 1 0 1 1 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 1 1 1 0 0 1 0 1 1 0 0 1 0 0 0 1 1 0 0 0 0 1 0 1 0 0 0 0 1 1 0 1 1 1 0 1 1 0 0 0 1 0 1 0 1 0 0 0 1 0 1 1 1 1 1 0 1 0 1 1 0 1 1 1 0 1 1 1 0 0 1 0 1 0 1 0 1 0 0 0 1 1 1 1 1 1 1 0 1 0 0 1 1 1 1 1 0 1 1 0 0 0 1 1 1 0 1 0 0 0 1 0 0 1 0 1 1 0 1 0 1 0 0 1 1 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 0 1 1 1 1 1 0 0 0 0 0 0 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0 1 1 0 0 0 1 0 1 1 0 0 0 0 1 0 1 1 0 0 0 1 1 0 1 1 0 0 0 1 0 0 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 0 1 0 1 0 1 0 1 1 0 1 1 1 1 1 1 1 0 0 1 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 0 1 0 0 0 1 0 1 0 0 0 1 0 1 0 1 1 1 0 1 0 1 1 1 1 0 1 0 1 1 1 0 1 1 0 1 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 1 0 0 1 1 1 1 0 1 1 0 0 0 0 1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 0 0 0 0 1 0 0 0 1 1 0 1 1 0 1 0 1 0 0 0 1 0 0 1 1 1 1 0 1 0 1 0 0 0 1 1 0 1 1 1 1 0 1 0 0 0 1 0 1 1 0 1 1 0 1 0 1 1 0 0 1 0 0 1 1 1 0 0 0 0 1 0 0 0 1 0 0 1 1 0 1 0 1 0 1 0 0 0 0 1 1 1 1 1 1 0 1 1 0 0 1 1 1 1 0 1 0 0 0 0 0 1 1 0 1 1 0 1 1 1 0 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 0 0 1 1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 1 0 1 1 0 1 0 0 1 1 1 1 0 0 1 1 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 0 1 0 1 1 0 0 1 1 0 1 1 1 0 0 0 0 0 0 1 0 1 0 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 1 1 1 0 0 1 0 1 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 1 0 0 1 0 0 1 1 1 0 1 0 0 1 0 0 0 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 0 1 1 1 1 0 0 0 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 1 0 1 0 0 1 0 1 1 1 0 1 1 0 0 1 1 0 1 0 1 0 0 0 0 0 0 1 1 1 1 0 1 1 1 1 0 0 1 1 0 1 0 1 1 0 0 0 0 0 1 1 1 0 0 0 1 1 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 1 0 1 0 0 1 0 0 1 1 1 1 1 0 0 1 0 0 0 1 1 1 0 0 0 1 0 1 0 0 1 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 1 0 0 0 1 0 0 0 0 1 1 0 1 0 1 0 1 1 0 0 0 1 1 1 1 1 0 0 0 1 0 0 1 1 1 0 0 1 0 1 0 0 0 1 0 0 0 1 1 1 0 1 0 1 0 1 0 0 1 0 1 1 1 1 1 1 0 0 1 1 0 1 1 1 1 0 0 0 0 0 1 0 1 1 0 0 1 1 1 0 1 1 0 0 0 0 0 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 0 0 1 0 1 0 1 0 1 1 0 0 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 1 1 0 1 0 1 1 1 1 0 0 0 1 1 1 0 1 1 0 0 1 0 0 1 0 1 0 0 0 0 1 0 0 1 1 1 0 1 1 0 1 0 0 0 1 0 1 0 0 1 1 0 1 0 1 1 1 0 0 0 0 1 1 1 0 1 0 0 1 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 0 0 0 1 1 1 1 0 0 0 0 1 0 0 1 1 0 0 1 1 0 1 0 0 0 0 0 0 0 0 1 1 0 1 1 1 1 1 1 0 0 0 1 0 1 1 1 1 0 0 1 0 1 1 0 1 1 0 0 0 1 1 0 0 1 0 0 0 1 0 0 0 0 0 1 0 1 0 1 0 1 1 1 0 1 1 1 1 1 1 0 1 0 1 0 1 1 1 1 0 1 1 1 1 1 0 1 1 0 1 0 1 1 1 0 1 0 0 1 1 1 0 1 0 1 1 0 0 0 1 0 1 1 1 0 0 0 1 0 1 1 0 1 0 0 1 0 1 1 0 0 1 1 0 0 1 1];
    
end

% ******** Sample the pattern based on the bit stream ********************

output_pattern = (output_pattern - 0.5)*2;

% Now we have the outptut bit pattern for a single PRBS.
% Construct the output signal.
patt = [];

for i = 1:length(output_pattern)
   patt = [patt, output_pattern(i)*ones(1,samples_per_bit)]; 
end

pattern = [];

for i = 1:num_repeats
    pattern = [pattern, patt];
end




end