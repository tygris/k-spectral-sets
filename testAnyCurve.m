%Any-curve branch test file

%Natalie Wellen
%10/25/21

%% Test findr 
A1 = [.5 -1; 1 -1.5];

sigma0 = -1/4 + 1/2*1i;
sigma0_prime = -1i;

[rA1, A1_1or2] = findr(A1, sigma0, sigma0_prime, 1, 1000)

%compare to previously tested plots
A1num_range = numerical_range(A1, 100);
spectral_set_choices(A1, A1num_range, 50);