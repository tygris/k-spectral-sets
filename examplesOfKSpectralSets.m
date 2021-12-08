%File of different examples of how K-Spectral Sets can be used

%Natalie Wellen
%12/07/21

%% Example 1
% Using the transient_demo() from Eigtool: https://www.cs.ox.ac.uk/pseudospectra/eigtool/

A =  transient_demo(10);

[nrA, nrA_prime] = numerical_range(A, 2000); 
figure()
plot(nrA), hold on
ewsA = eig(A);
plot(ewsA, 'kx')

pause

% The numerical range  can be reduced to bound the transients
rhm = -10^-8; %right-hand max
del_OmA = nrA;
del_OmA_prime = -1i*(del_OmA +1/2);
vertical_index = (real(del_OmA)> rhm);
del_OmA(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_OmA(vertical_index));
del_OmA_prime(vertical_index) = -1i;

%calculate spectral set value
[kA, c1A, c2A] = calc_k(A, nrA, nrA_prime, del_OmA, del_OmA_prime, 1, 1, 1000)
close

% calculate matrix envelope
timestep = 0.1; iterations = 500;
meA = times_expm(A, timestep, iterations);
t = 0:timestep:iterations*timestep;

figure()
plot(t, kA*exp(t*(-10^-8)), 'DisplayName', 'K'), hold on
plot(t, meA)

pause

%Now we compare the value of the K to the corresponding pseudospectral
%   bounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First very small dimension systems

% for jj = 3:9
%     A = transient_demo(jj);
%     eigtool(A)
% end
%all of the relevant information is then exported and saved to td_39.mat
load("td_39.mat")
%%
dims = 4:9; dimlength = 6;
Gam_tds = cell(1,dimlength); td_lbs = cell(1,dimlength); td_ups = cell(1,dimlength);

for jj = 1:dimlength
    Gam_td = pe_contour(cell2mat(xtds(jj)),cell2mat(ytds(jj)),cell2mat(Ztds(jj)),cell2mat(powertds(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :))
    Gam_tds(jj) = {Gam_td};
    td_lbs(jj) = {td_lb};
    td_ups(jj) = {td_up};
end

% Now find our upper bound on the matrix envelope
td_ks = cell(1,dimlength); td_c1s = cell(1,dimlength); td_c2s = cell(1,dimlength);
for jj = 1:dimlength
    A = transient_demo(jj+3);
    [nr, nr_prime] = numerical_range(A, 500);
    rhm = -10^-8; %right-hand max
    del_Om = nr;
    del_Om_prime = -1i*(del_Om +1/2);
    vertical_index = (real(del_Om)> rhm);
    del_Om(vertical_index) = rhm*ones(1,sum(vertical_index)) + 1i*imag(del_Om(vertical_index));
    del_Om_prime(vertical_index) = -1i;

    %calculate spectral set
    [k, c1, c2] = calc_k(A, nr, nr_prime, del_Om, del_Om_prime, 1, 1, 500) 
    td_ks(jj) = {k};
    td_c1s(jj) = {c1};
    td_c2s(jj) = {c2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Then larger dense systems
load("transientdem_dimex.mat")
Gam_tdl = cell(1,11); td_lbl = cell(1,11); td_upl = cell(1,11);
dims = 10:4:50;
for jj = 1:6
    Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :))
    Gam_tdl(jj) = {Gam_td};
    td_lbl(jj) = {td_lb};
    td_upl(jj) = {td_up};
end

% Now find our upper bound on the matrix envelope
td_kl = cell(1,11); td_c1l = cell(1,11); td_c2l = cell(1,11);
for jj = 1:6
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
    td_kl(jj) = {k};
    td_c1l(jj) = {c1};
    td_c2l(jj) = {c2};
end
    
    
    
    
    
    
    
    
    
    
    