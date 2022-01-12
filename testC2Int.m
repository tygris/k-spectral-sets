%Script with tests for c2-int

%Natalie Wellen
%1/06/22

%% Test calc_c2 
% 1. With A1 and the right half plane removed
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
figure(), plot(del_OmA1), daspect([1,1,1])

%[c2A1] = calc_c2(A1, nrA1, del_OmA1, del_OmA1_prime, 1, 100)
[kA1, c1A1, c2A1] = calc_k(A1, nrA1, del_OmA1, del_OmA1_prime, 1, 100)

pause 

% 2. Test calc_c2 with A1 and two regions removed that exclude nrA1(1)
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
vertical_index = (real(del_OmA1)< -1);
del_OmA1(vertical_index) = -1*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = 1i;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 100)

pause

% 3. Test calc_c2 with A1 and two regions removed that don't exclude nrA1(1)
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
horizontal_index = (imag(del_OmA1)> 0.75);
del_OmA1(horizontal_index) = 0.75*ones(1,sum(horizontal_index))*1i + real(del_OmA1(horizontal_index));
del_OmA1_prime(horizontal_index) = 1;
horizontal_index = (imag(del_OmA1)< -0.75);
del_OmA1(horizontal_index) = -0.75*ones(1,sum(horizontal_index))*1i + real(del_OmA1(horizontal_index));
del_OmA1_prime(horizontal_index) = -1;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 100)







%% Compare calc_c2_old and calc_c2 results

%First using the demo I created
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 2000);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -10^-8);
del_OmA1(vertical_index) = -10^-8*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
figure(), plot(del_OmA1), daspect([1,1,1])

tic, c2_new = calc_c2(A1, nrA1, del_OmA1, del_OmA1_prime, 1, 500), toc
tic, c2_old = calc_c2_old(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 500), toc

% c2_new took about 23 minutes at 1354.583707 seconds for a result of 1.7684
% c2_old took 5.133408 seconds for a result of 2.7793
% k_new = 3.4322 and k_old = 4.7233
%Conclusions: I was using this example specifically because it was well
% behaved, and I think that played a role in the small gains for c2_new.
% It may also be the case that c2_new would need much fewer points and 
% coud run similarly fast to get the same result/error. 
% I also think that if I chose to compare on the rainforest matrix or
% the Boeing demo matrix from eigtool, then I would find c2_new worth the
% wait. I test this hypothesis next with the rainforest matrix.


%% Second using the rainforest matrix
A =  [-1.5622, 0.6685,0,0,0,0,0,0,0;
     0, -0.7119,0,0,2.5632,0,0,0,0;
     1.4627, 0.0364, -6.4091,0,0,1.1446,0,55.8201,17.2972;
     0,0,0,-0.0222,0,0,315.9443,0,0;
     0,0,0,0.0201,-2.5632,0,0,0,0;
     0,0.0070,0,0,0,-2.0348,0,0,0;
     0,0,6.4091,0,0,0,-315.9443,0,0;
     0.0995,0,0,0,0,0.8902,0,-62.6458,0;
     0,0,0,0,0,0,0,6.8257,-17.2972];

[nrA, nrA_prime] = numerical_range(A, 400); 
% figure()
% plot(nrA), hold on
% ewsA = eig(A);
% plot(ewsA, 'kx')

% The numerical range  can be reduced to bound the transients
rhm = 0; %right-hand max
del_OmA = nrA;
del_OmA_prime = nrA_prime;
vertical_index = (real(del_OmA)> rhm);
del_OmA(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA(vertical_index));
del_OmA_prime(vertical_index) = -1i;

tic, calc_c2_old(A, nrA, nrA_prime, del_OmA, del_OmA_prime, 1, 1000), toc
tic, calc_c2(A, nrA, del_OmA, del_OmA_prime, 1, 1000), toc

% calc_c2 took about 7126.5 seconds to run and got a value of
% 350.2867730261789
% calc_c2_old took 36.27 seconds to run, and got a value of
% 56866.7472436604
% As expected the benefits of the new method are much more noticeable here.
% We were able to reduce the upper bound on c2 by over 162 times.
% However, the result of c2<=350.3 is nowhere near the value of 7 I was expecting based on
% conversations with my advisor.


%% Look for order of convergence of the method
% I will use the transient demo matrix from eigtool
% First by increasing the number of points defining the boundary
A = transient_demo(10);
rhm = 0; %right-hand max

A10_res100_c2 = zeros(3,6);
count = 1;
for jj = [100, 200, 400, 800, 1600, 3200]
    [nrA, nrA_prime] = numerical_range(A, jj);
    del_OmA = nrA;
    del_OmA_prime = nrA_prime;
    vertical_index = (real(del_OmA)> rhm);
    del_OmA(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA(vertical_index));
    del_OmA_prime(vertical_index) = -1i;
    
    jj
    tic, c2 = calc_c2(A, nrA, del_OmA, del_OmA_prime, 1, 100), toc
    A10_res100_c2(1:2,count) = [jj; c2];
    count = count+1;
end


%% Second by increasing the resolution of the findr() function

A10_nr100_c2 = zeros(3,4);
[nrA, nrA_prime] = numerical_range(A, 100);
del_OmA = nrA;
del_OmA_prime = -1i*(del_OmA +1/2);
vertical_index = (real(del_OmA)> rhm);
del_OmA(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA(vertical_index));
del_OmA_prime(vertical_index) = -1i;

for jj = 1:4
    10^jj
    tic, c2 = calc_c2(A, nrA, del_OmA, del_OmA_prime, 1, 10^jj), toc
    A10_nr100_c2(1:2,jj) = [10^jj; c2];
end

%% Compare to calculating the resolvent norm along the boundary
% Question, is the resolvent norm bounded on the boundary of W(A)?

%My toy example
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define Gam1
rhm = 0;
vertical_index = (real(nrA1)> rhm);
Gam1 = 1i*imag(nrA1(vertical_index));
Gam1_prime(vertical_index) = -1i;
%calculate the cif on Gam1
cif = cauchyIntFormula(A1, Gam1);

%calculate the cif for the boundary of Omega
delOmA1 = nrA1; 
delOmA1(vertical_index) = 1i*imag(delOmA1(vertical_index));
cifA1 = cauchyIntFormula(A1, delOmA1)

pause
%
% calculate matrix envelope
timestep = 0.1; iterations = 500;
meA = times_expm(A1, timestep, iterations);
t = 0:timestep:iterations*timestep;
K = 3.432250954863446; %using 1.7684 calc from first section for c2

figure()
plot(t, cifA1*ones(1,iterations+1), 'DisplayName', 'CIF'), hold on
plot(t, K*ones(1, iterations+1), 'DisplayName', 'K')
plot(t, meA, 'DisplayName', 'exp(At)')
legend()

%% Compare to calculating the resolvent norm along the boundary
% The Panamanian Rainforest Example

A2 =  [-1.5622, 0.6685,0,0,0,0,0,0,0;
     0, -0.7119,0,0,2.5632,0,0,0,0;
     1.4627, 0.0364, -6.4091,0,0,1.1446,0,55.8201,17.2972;
     0,0,0,-0.0222,0,0,315.9443,0,0;
     0,0,0,0.0201,-2.5632,0,0,0,0;
     0,0.0070,0,0,0,-2.0348,0,0,0;
     0,0,6.4091,0,0,0,-315.9443,0,0;
     0.0995,0,0,0,0,0.8902,0,-62.6458,0;
     0,0,0,0,0,0,0,6.8257,-17.2972];

[nrA2] = numerical_range(A, 400000); 
% The numerical range  can be reduced to bound the transients
rhm = 0; %right-hand max
del_OmA2 = nrA2;
vertical_index = (real(del_OmA2)> rhm);
del_OmA2(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA2(vertical_index));

cifA3 = cauchyIntFormula(A2, del_OmA2)

pause

% calculate matrix envelope
timestep = 0.1; iterations = 500;
meA2 = times_expm(A2, timestep, iterations);
t = 0:timestep:iterations*timestep;
K2 = 350.2867730261789 + sqrt(351.2867730261789); %using calc from second section for c2

figure()
plot(t, cifA3*ones(1,iterations+1), 'DisplayName', 'CIF'), hold on
plot(t, K2*ones(1, iterations+1), 'DisplayName', 'K')
plot(t, meA2, 'DisplayName', 'exp(At)')
legend()

pause
% Calculate the resolvent norm just for A2's Gam1
breaks = find(vertical_index(2:end) - vertical_index(1:end-1));
Gam1 = cat(2,del_OmA2(breaks(2)+1:end), del_OmA2(1:breaks(1)));

cif = cauchyIntFormula(A3, Gam1)



%% Compare to calculating the resolvent norm along the boundary
% The transient demo of dimension 10 from the eigtool package

A3 = transient_demo(10);
[nrA3] = numerical_range(A, 400); 
% The numerical range  can be reduced to bound the transients
rhm = 0; %right-hand max
del_OmA3 = nrA3;
vertical_index = (real(del_OmA3)> rhm);
del_OmA3(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA3(vertical_index));

cifA3 = cauchyIntFormula(A3, del_OmA3)

pause

% calculate matrix envelope
timestep = 0.1; iterations = 500;
meA3 = times_expm(A3, timestep, iterations);
t = 0:timestep:iterations*timestep;
K3 = 3.0799+ sqrt(4.0799); %using calc from earlier section for c2

figure()
plot(t, cifA3*ones(1,iterations+1), 'DisplayName', 'CIF'), hold on
plot(t, K3*ones(1, iterations+1), 'DisplayName', 'K')
plot(t, meA3, 'DisplayName', 'exp(At)')
legend()

pause 
% Calculate the resolvent norm just for A3's Gam1
breaks = find(vertical_index(2:end) - vertical_index(1:end-1));
Gam1 = cat(2,del_OmA3(breaks(2)+1:end), del_OmA3(1:breaks(1)));

cif = cauchyIntFormula(A3, Gam1)



































