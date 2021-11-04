%Any-curve branch test file

%Natalie Wellen
%10/25/21

%% Test findr 
A1 = [.5 -1; 1 -1.5];

sigma0 = -1/4 + 1/2*1i;
sigma0_prime = -1i;

[rA1, A1_1or2] = findr(A1, sigma0, sigma0_prime, 1, 1000)

%compare to previously tested plots
[nrA1, nrA1_prime] = numerical_range(A1, 100);
%spectral_set_choices(A1, nrA1, 50);

%% Test calc_c2 with A1
% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1, minr_A1, A1_1or2] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 1000)
 %for this example the minimum value of c2 happens on the real line (del_OmA1(1))
 
 % Calculate c1 for A1 and this reshaped boundary
 c1A1 = calc_c1(1, -pi/2, del_OmA1)
 c1A1 = calc_c1(1, 3*pi/2, del_OmA1)
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
 
% calculate k for A1 and this reshaped boundary 
 
%k = calc_k(A1, del_OmA1, del_OmA1_prime, 1, 1, 1000)

%% What happens when the numerical range is passed for calculating c2 (nothing removed)?

[c2A1_2, minr_A1_2, A1_1or2_2] = calc_c2(A1, nrA1, nrA1_prime, nrA1, -1i*(del_OmA1 +1/2), 1, 1000)

 %% Find c2 for the complicated Panamanian Rainforest example from Neubert and Caswell
 
A2 = [-1.5622, 0.6685,0,0,0,0,0,0,0;
     0, -0.7119,0,0,2.5632,0,0,0,0;
     1.4627, 0.0364, -6.4091,0,0,1.1446,0,55.8201,17.2972;
     0,0,0,-0.0222,0,0,315.9443,0,0;
     0,0,0,0.0201,-2.5632,0,0,0,0;
     0,0.0070,0,0,0,-2.0348,0,0,0;
     0,0,6.4091,0,0,0,-315.9443,0,0;
     0.0995,0,0,0,0,0.8902,0,-62.6458,0;
     0,0,0,0,0,0,0,6.8257,-17.2972];
A2_eigs = eig(A2);

A2_fov = fov(A2);
A2_fov_prime = diff(A2_fov);

% convert to discrete vectors with nump points
L = 2*pi;
nump = 1000;
ds = L/(nump);
nrA2 = A2_fov([0:nump]*ds);
nrA2_prime = A2_fov_prime([0:nump]*ds);
nrA2 = flip(nrA2);
nrA2_prime = flip(nrA2_prime);

% reshape del_Omega so that the real part is always less than zero
delta = 1e-6;
del_OmA2 = nrA2; del_OmA2_prime = nrA2_prime;
vertical_index = (real(del_OmA2)> -delta);
del_OmA2(vertical_index) = -delta*ones(1,sum(vertical_index)) + 1i*imag(del_OmA2(vertical_index));
del_OmA2_prime(vertical_index) = -1i;
figure(), plot(del_OmA2), daspect([1,1,1])
hold on
plot(A2_eigs, zeros(1, 9),'x')

% Calculate c1 and the min eigenvalues for A2

c1A2 = calc_c1(1, -pi/2, del_OmA2)
[c2A2, minr_A2, A2_1or2] = calc_c2(A2, nrA2, nrA2_prime, del_OmA2, del_OmA2_prime, .1, 1000)
 
 
%% Test curve creation for measuring arclength in calc_c2

A1 = [.5 -1; 1 -1.5];
nrA1 = numerical_range(A1, 100);
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
fova = nrA1;
del_Om = del_OmA1;
in1 = ~ismember(fova, del_Om);
in2 = ~ismember(del_Om, fova);
out = fova(in1);
in = fova(~in1);
figure(), plot(out), daspect([1,1,1]), hold on
plot(out(1:4), 'o')
plot(in(1), 'o')

%%
in1 = ~ismember(nr, del_Om);
in2 = ~ismember(del_Om, nr);
if in1(1)
    ind1 = in1 & imag(nr)>=0;
    ind2 = in2 & imag(del_Om)>=0;
    ind3 = in1 & imag(nr) < 0;
    ind4 = in2 & imag(del_Om) < 0;
    gam1 = cat(2, nr(ind1), flip(del_Om(ind2)), flip(del_Om(ind4)), nr(ind3), nr(1));
%    gam1_prime = cat(2, nr_prime(ind1), flip(del_Om_prime(ind2)), flip(del_Om_prime(ind4)), nr_prime(ind3), nr_prime(1));
    figure(), plot(gam1), hold on
else
    temp = nr(in1);
    tempp = nr_prime(in1);
    gam1 = cat(2, temp, flip(del_Om(in2)), temp(1));
%    gam1_prime = cat(2, tempp, flip(del_Om_prime(in2)), tempp(1));
    figure(), plot(gam1), hold on
end
%calculate the arclength of gam1 using the absolute distance between
%points
dz = abs(gam1(2:end) - gam1(1:end-1));
L = sum(dz)

 
 
 
 
 