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

%% Test calc_c2 with A1
% define del_Om for testing
del_OmA1 = A1num_range;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1, minr_A1, A1_1or2] = calc_c2(A1, del_OmA1, del_OmA1_prime, 1, 1000)
 %for this example the minimum value of c2 happens on the real line (del_OmA1(1))
 
 %% Calculate c1 for A1 and this reshaped boundary
 c1A1 = c1_estimate(1, -pi/2, del_OmA1)
 c1A1 = c1_estimate(1, 3*pi/2, del_OmA1)
 %UHOH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 %first error is symmetric across the break point gets a zero:
 a1 = angle(del_OmA1(50) - del_OmA1(1)); 
 a2 = angle(del_OmA1(51) - del_OmA1(1));
 theta = min(abs(a1+a2), abs(a1 - a2));
 
 %second error is that a positive angle in the negative imaginary axis
 %is assumed to be passed for sigma0_prime
 %angle_0 in c1_estimate
 angle_0 = mod(angle(del_OmA1(2)-del_OmA1(1)), 2*pi);
 angle_0 = min(abs(angle_0 - mod(-pi/2+pi, 2*pi)), abs(-pi/2 -angle_0));
 
%% calculate k for A1 and this reshaped boundary 
 
k = calc_k(A1, del_OmA1, del_OmA1_prime, 1, 1, 1000)

%% What happens when the numerical range is passed for calculating c2 (nothing removed)?

[c2A1_2, minr_A1_2, A1_1or2_2] = calc_c2(A1, A1num_range, -1i*(del_OmA1 +1/2), 1, 1000)

 %% Find c2 for the complicated Panamanian Rainforest example from Neubert and Caswell
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 