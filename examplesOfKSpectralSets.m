%File of different examples of how K-Spectral Sets can be used

%Natalie Wellen
%12/07/21

%% Example 1
% Using the transient_demo() from Eigtool: https://www.cs.ox.ac.uk/pseudospectra/eigtool/

A =  transient_demo(10);

[k, cif] = contDS(A,2)


%% Now we compare the value of the K to the corresponding pseudospectral
%   bounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First very small dimension systems

% for jj = 3:9
%     A = transient_demo(jj);
%     eigtool(A)
% end
%all of the relevant information is then exported and saved to td_39.mat
load("td_39.mat")
%
dims = 4:9; dimlength = 6;
Gam_tds = cell(1,dimlength); td_lbs = cell(1,dimlength); td_ups = zeros(1,dimlength);

for jj = 1:dimlength
    Gam_td = pe_contour(cell2mat(xtds(jj)),cell2mat(ytds(jj)),cell2mat(Ztds(jj)),cell2mat(powertds(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c')
    td_up = exp(1)*dims(jj)*max(td_lb(2, :))
    Gam_tds(jj) = {Gam_td};
    td_lbs(jj) = {td_lb};
    td_ups(jj) = td_up;
end

% Now find our upper bound on the matrix envelope
td_ks = zeros(1,dimlength); td_cifs = zeros(1,dimlength); 
td_c1s = zeros(1,dimlength); td_c2s = zeros(1,dimlength);
for jj = 1:dimlength
    A = transient_demo(jj+3);
    figure()
    [k,cif, c1,c2] = contDS(A,2);
    td_ks(jj) = k;
    td_cifs(jj) = cif;
    td_c1s(jj) = c1;
    td_c2s(jj) = c2;
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
    [k, c1, c2] = calc_k(A, nr, nr_prime, del_Om, del_Om_prime, 1, 500) 
    td_kl(jj) = {k};
    td_c1l(jj) = {c1};
    td_c2l(jj) = {c2};
end
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Using the contDS Function pass these matrices

A = [.5 -1; 1, -1.5];
figure()
[kA, cifA] = contDS(A)

C = [0 1 2; -0.01 0 3; 0 0 0];
figure()
[kC, cifC] = contDS(C)

T = [-0.9503,0,0.0690,0.0002, 0.0027,0.0034;
         0.95, -0.18,0,0,0,0;
         0, 0.15, -0.2569,0,0,0;
         0,0,0.100, -0.0138,0,0;
         0,0,0.0019, 0.0002, -0.0124,0;
         0,0,0,0.0001, 0.0028, -0.0049]; %1986 Tuesday Lake
figure()
[kT,cifT] = contDS(T)

RF = zeros(9);
RF(1,1) = -1.5622; RF(1,2) = .6685;
RF(2,2) = -.7119; RF(2,5) = 2.5632;
RF(3,1) = 1.4627; RF(3,2) = .0364; RF(3,3) = -6.4091; RF(3,6) = 1.1446; RF(3,8) = 55.8201;
RF(3,9) = 17.2972;
RF(4,4) = -.0222; RF(4,7) = 315.9443;
RF(5,4) = .0201; RF(5,5) = -2.5632;
RF(6,2) = .0070; RF(6,6) = -2.0348;
RF(7,3) = 6.4091; RF(7,7) = -315.9443;
RF(8,1) = .0995; RF(8,6) = .8902; RF(8,8) = -62.6458;
RF(9,8) = 6.8257; RF(9,9) = -17.2972;
figure()
[kRF, cifRF] = contDS(RF)

B = boeing_demo('S');
figure()
[kB, cifB] = contDS(B,2);

   




































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block Diagonal Matrix Example
    
b1 = 3;
b2 = 5;
b3 = 2;
res = 10000;
A = [randn(b1), zeros(b1, b2); zeros(b2, b1), randn(b2)+(10+5i)*eye(b2)];%, zeros(b2, b3); zeros(b3, b1+b2), rand(b3)-10-10*1i];

figure()
nr = numerical_range(A, res);
plot(nr), hold on
evs = eig(A);
plot(evs, 'kx')
plot(numerical_range(A(1:b1,1:b1), res))
plot(numerical_range(A(b1+1:b1+b2,b1+1:b1+b2), res))
%plot(numerical_range(A(b1+b2+1:b1+b2+b3,b1+b2+1:b1+b2+b3), res))



