%Any-curve branch test file

%Natalie Wellen
%12/01/21

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

[c2A1, mineigA1, LA1] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 1000)
 %for this example the minimum value of c2 happens on the real line (del_OmA1(1))
 
 %% compare the calculated arclength to the exact arclength
 % to do this take advantage of the nr of A1 being a circle with radius 1
 index = ismember(del_OmA1, nrA1);
 temp = del_OmA1(index);
 l1 = angle_between(temp(1) + 0.5, temp(end)+0.5);
 l2 = abs(temp(1) - temp(end));
 L_compare = l1 + l2;
 accuracy_arclength = abs(LA1-L_compare)
 relative_accuracy_arclength = accuracy_arclength/LA1
 
 %% compare the calculated arclength to the exact arclength
 % to do this take advantage of the nr of A1 being a circle with radius 1
 index = ismember(del_OmA1, nrA1);
 temp = del_OmA1(~index);
 temp1 = max(imag(temp)); temp2 = min(imag(temp));
 l1 = angle_between(0.25+1i*temp1, 1i*temp2+0.25);
 l2 = abs(temp1 - temp2);
 L_compare = l1 + l2;
 accuracy_arclength = abs(LA1-L_compare)
 relative_accuracy_arclength = accuracy_arclength/LA1
 
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

%% Check accuracy of measureArcLength

% Start with measuring the arcLength of a circle
accuracy = zeros(1,7);
relative_accuracy = zeros(1,7);
for ii = 1:7
    circ = circle(1,0,10^ii);
    L_temp = measureArcLength(circ, NaN);
    accuracy(ii) = abs(2*pi-L_temp);
    relative_accuracy(ii) = accuracy(ii)/2*pi;
end
accuracy, relative_accuracy

%% What happens when the numerical range is passed for calculating c2 (nothing removed)?

[c2A1_2, minr_A1_2, A1_1or2_2] = calc_c2(A1, nrA1, nrA1_prime, nrA1, -1i*(del_OmA1 +1/2), 1, 1000)
























%% Example 1
A1 = [.5 -1; 1 -1.5];

% Using matrix A1 from before but with higher resolution.
[nrA1, nrA1_prime] = numerical_range(A1, 2000); 
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -10^-8);
del_OmA1(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;

%calculate spectral set
[kA1, c1A1, c2A1] = calc_k(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 1, 1000)
close

%calculate matrix envelope
timestep = 0.01; iterations = 500;
meA1 = times_expm(A1, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA1*exp(t*(-10^-8)), 'DisplayName', 'K')
hold on, daspect([1,1,1])
plot(t, meA1)

%% Example 2
 
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

k = c2A2 + sqrt(c1A2^2 + c2A2);                                                                                                                                                                                                                                                                                                                                                       
spectral_upper_bound = exp(-delta)*k
pseudospectral_upper_bound = 64.803447271835608 %see CaswellNeubert_PseudoExp.m in the pseudospectral_bounds package for calculation.
 


%% Example 3; Interesting because I can't calculate it, but eigtool can...
A3 = boeing_demo('s');

[nrA3, nrA3_prime] = numerical_range(A3, 2000);
%
del_OmA3 = nrA3;
del_OmA3_prime = -1i*(del_OmA3 +1/2);
vertical_index = (real(del_OmA3)> -10^-8);
del_OmA3(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA3(vertical_index));
del_OmA3_prime(vertical_index) = -1i;

%calculate spectral set
[kA3, c1A3, c2A3] = calc_k(A3, nrA3, nrA3_prime, del_OmA3, del_OmA3_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.01; iterations = 2000;
meA3 = times_expm(A3, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
%plot(t, kA3*exp(t*(-10^-8)), 'DisplayName', 'K')
%hold on, daspect([1,1,1])
plot(t, meA3)



%% Example 4

A4 = markov_demo(50);

[nrA4, nrA4_prime] = numerical_range(A4, 2000); 
figure()
plot(nrA4), hold on
ewsA4 = eig(A4);
plot(ewsA4, 'kx')

%
del_OmA4 = nrA4;
del_OmA4_prime = -1i*(del_OmA4 +1/2);
vertical_index = (real(del_OmA4)> -10^-8);
del_OmA4(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA4(vertical_index));
del_OmA4_prime(vertical_index) = -1i;

%calculate spectral set
[kA4, c1A4, c2A4] = calc_k(A4, nrA4, nrA4_prime, del_OmA4, del_OmA4_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.01; iterations = 500;
meA4 = times_expm(A4, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA4*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA4)



%% %% Example 5
A5 =  airy_demo(10);

[nrA5, nrA5_prime] = numerical_range(A5, 2000); 
figure()
plot(nrA5), hold on
ewsA5 = eig(A5);
plot(ewsA5, 'kx')
% This numerical range is already entirely in the left-half plane
del_OmA5 = nrA5;
del_OmA5_prime = -1i*(del_OmA5 +1/2);
vertical_index = (real(del_OmA5)> -10^-8);
del_OmA5(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA5(vertical_index));
del_OmA5_prime(vertical_index) = -1i;

%calculate spectral set
[kA5, c1A5, c2A5] = calc_k(A5, nrA5, nrA5_prime, del_OmA5, del_OmA5_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.01; iterations = 500;
meA5 = times_expm(A5, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA5*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA5)



%% %% Example 6
A6 =  basor_demo(10);

[nrA6, nrA6_prime] = numerical_range(A6, 2000); 
figure()
plot(nrA6), hold on
ewsA6 = eig(A6);
plot(ewsA6, 'kx')
% this has a spectral set that cannot be contained in the left-half plane
% nor the unit disk
del_OmA6 = nrA6;
del_OmA6_prime = -1i*(del_OmA6 +1/2);
vertical_index = (real(del_OmA6)> -10^-8);
del_OmA6(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA6(vertical_index));
del_OmA6_prime(vertical_index) = -1i;

%calculate spectral set
[kA6, c1A6, c2A6] = calc_k(A6, nrA6, nrA6_prime, del_OmA6, del_OmA6_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.01; iterations = 500;
meA6 = times_expm(A6, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA6*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA6)



%% %% Example 7
A7 =  chebspec_demo(10);

[nrA7, nrA7_prime] = numerical_range(A7, 2000); 
figure()
plot(nrA7), hold on
ewsA7 = eig(A7);
plot(ewsA7, 'kx')
% This numerical range can be dramatically reduced... but doesn't fit the ds application
rhm = 10; %right-hand max
del_OmA7 = nrA7;
del_OmA7_prime = -1i*(del_OmA7 +1/2);
vertical_index = (real(del_OmA7)> rhm);
del_OmA7(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA7(vertical_index));
del_OmA7_prime(vertical_index) = -1i;

%calculate spectral set
[kA7, c1A7, c2A7] = calc_k(A7, nrA7, nrA7_prime, del_OmA7, del_OmA7_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.01; iterations = 500;
meA7 = times_expm(A7, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA7*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA7)

%from our work, 118.35*f(10+ ix) x in [22.6i, -22,6i]
%general numerical range says 1+ sqrt(2) f(32)


%% %% Example 8
A8 =  transient_demo(10);

[nrA8, nrA8_prime] = numerical_range(A8, 2000); 
figure()
plot(nrA8), hold on
ewsA8 = eig(A8);
plot(ewsA8, 'kx')
%% This numerical range can be dramatically reduced... but doesn't fit the ds application
rhm = -10^-8; %right-hand max
del_OmA8 = nrA8;
del_OmA8_prime = -1i*(del_OmA8 +1/2);
vertical_index = (real(del_OmA8)> rhm);
del_OmA8(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA8(vertical_index));
del_OmA8_prime(vertical_index) = -1i;

%calculate spectral set
[kA8, c1A8, c2A8] = calc_k(A8, nrA8, nrA8_prime, del_OmA8, del_OmA8_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.1; iterations = 500;
meA8 = times_expm(A8, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA8*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA8)



%% compare pseudospectral upper bounds and ours as the dimension of the matrix increases
% I'm not seeing anything that impressive in small dimensions, maybe the
% "gold" is in high dimensional systems?

% for jj = [10:4:50]
%     A = transient_demo(jj);
%     eigtool(A);
% end
load("transientdem_dimex.mat")

Gam_tds = cell(1,11); td_lbs = cell(1,11); td_ups = cell(1,11);
dims = 10:4:50;
for jj = 1:11
    Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :))
    Gam_tds(jj) = {Gam_td};
    td_lbs(jj) = {td_lb};
    td_ups(jj) = {td_up};
end

%% Now find our upper bound on the matrix envelope
td_k = cell(1,11); td_c1 = cell(1,11); td_c2 = cell(1,11);
for jj = 1:11
    A = transient_demo(4*jj+6);
    [nr, nr_prime] = numerical_range(A, 500);
    rhm = -10^-8; %right-hand max
    del_Om = nr;
    del_Om_prime = -1i*(del_Om +1/2);
    vertical_index = (real(del_Om)> rhm);
    del_Om(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_Om(vertical_index));
    del_Om_prime(vertical_index) = -1i;

    %calculate spectral set
    [k, c1, c2] = calc_k(A, nr, nr_prime, del_Om, del_Om_prime, 1, 1, 500) 
    td_k(jj) = {k};
    td_c1(jj) = {c1};
    td_c2(jj) = {c2};
end








                                   













 
 
 
 
 