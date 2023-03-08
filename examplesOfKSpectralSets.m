%File of different examples of how K-Spectral Sets can be used

%Natalie Wellen
%2/10/23

%% Example 1
% Using the transient_demo() from Eigtool: https://www.cs.ox.ac.uk/pseudospectra/eigtool/

A =  transient_demo(10);

[k, cif, c2] = contDS(A,2)

[k, cif, c2] = discDS(A, .75) %ylim([1, 12])

%% Explore different sets for calculating the integral of the resolvent norm

cif_i = resolvent_norm_integral(A, 1i*(-10000:0.01:10000))

DU = @(x,r) r*exp(1i*x);
unitD = DU(0:0.01:2*pi,1);
cif_ud = resolvent_norm_integral(A, unitD)

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
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :));
    Gam_tds(jj) = {Gam_td};
    td_lbs(jj) = {td_lb};
    td_ups(jj) = td_up;
end

% Now find our upper bound on the matrix envelope
td_ks = zeros(1,dimlength); td_cifs = zeros(1,dimlength); 
td_c2s = zeros(1,dimlength);
for jj = 1:dimlength
    A = transient_demo(jj+3);
    [k,cif,c2] = contDS(A,2);
    td_ks(jj) = k;
    td_cifs(jj) = cif;
    td_c2s(jj) = c2;
end

td_ups
td_ks
td_cifs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Then larger dense systems
 load("transientdem_dimex.mat")
 dims = 10:4:50; 
 dimlength = length(dims);
Gam_tdl = cell(1,dimlength); td_lbl = cell(1,dimlength); td_upl = zeros(1,dimlength);
for jj = 1:dimlength
    Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :));
    Gam_tdl(jj) = {Gam_td};
    td_lbl(jj) = {td_lb};
    td_upl(jj) = td_up;
end
%38 should be K(A) = 86.8702, eps = 10^-3.3
% Now find our upper bound on the matrix envelope
td_kl = zeros(1,dimlength); td_cifl = zeros(1,dimlength); td_c2l = zeros(1,dimlength);
for jj = 1:dimlength
    A = transient_demo(4*jj+6);
    %calculate spectral set
    [k, cif, c2] = contDS(A,4); 
    td_kl(jj) = k;
    td_cifl(jj) = cif;
    td_c2l(jj) = c2;
end
    
td_upl
td_kl
td_cifl

figure()
dims = [4:9, 10:4:50];
opts = {'LineWidth', 2};
semilogy(dims, [td_ups td_upl], 'DisplayName', 'Kreiss', opts{:})
hold on
plot(dims, [td_ks td_kl], 'DisplayName', 'K', opts{:})
plot(dims, [td_cifs td_cifl], 'DisplayName', 'Resolvent Norm', opts{:})
title('Transient Matrix Growth Bounds')
xlabel('dim(A)')
ylabel('Upper Bound of ||f(A)||')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transient Demo Discrete Version 
%          
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
    td_lb = pseudo_lb(Gam_td, 'd');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :));
    Gam_tds(jj) = {Gam_td};
    td_lbs(jj) = {td_lb};
    td_ups(jj) = td_up;
end

% Now find our upper bound on the matrix envelope
td_ks = zeros(1,dimlength); td_cifs = zeros(1,dimlength); 
td_c2s = zeros(1,dimlength);
for jj = 1:dimlength
    A = transient_demo(jj+3);
    [k,cif,c2] = discDS(A,2);
    td_ks(jj) = k;
    td_cifs(jj) = cif;
    td_c2s(jj) = c2;
end

td_ups
td_ks
td_cifs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Then larger dense systems
% load("transientdem_dimex.mat")
% dims = 10:4:50; dimlength = length(dims);
% Gam_tdl = cell(1,dimlength); td_lbl = cell(1,dimlength); td_upl = zeros(1,dimlength);
% for jj = 1:dimlength
%     Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
%     td_lb = pseudo_lb(Gam_td, 'd');
%     td_up = exp(1)*dims(jj)*max(td_lb(2, :));
%     Gam_tdl(jj) = {Gam_td};
%     td_lbl(jj) = {td_lb};
%     td_upl(jj) = td_up;
% end
% 
% % Now find our upper bound on the matrix envelope
% td_kl = zeros(1,dimlength); td_cifl = zeros(1,dimlength); td_c2l = zeros(1,dimlength);
% for jj = 1:dimlength
%     A = transient_demo(4*jj+6);
%     %calculate spectral set
%     [k, cif, c2] = discDS(A,4); 
%     td_kl(jj) = k;
%     td_cifl(jj) = cif;
%     td_c2l(jj) = c2;
% end
load("discreted_l.mat")    
td_upl
td_kl
td_cifl

figure()
dims = [4:9, 10:4:50];
opts = {'LineWidth', 2};
semilogy(dims, [td_ups td_upl], 'DisplayName', 'Kreiss', opts{:})
hold on
plot(dims, [td_ks td_kl], 'DisplayName', 'K', opts{:})
plot(dims, [td_cifs td_cifl], 'DisplayName', 'Resolvent Norm', opts{:})
%title('Transient Matrix Growth Bounds')
xlabel('dim(A)')
ylabel('Upper Bound of ||f(A)||')
legend()
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Using the contDS Function, pass these matrices

A = [.5 -1; 1, -1.5];
[kA, cifA, c2, cifG] = contDS(A)

C = [0 1 2; -0.01 0 3; 0 0 0];
[kC, cifC] = discDS(C)
%C is a bad example for what we do because one of its eigenvalues = 0+0i

T = [-0.9503,0,0.0690,0.0002, 0.0027,0.0034;
         0.95, -0.18,0,0,0,0;
         0, 0.15, -0.2569,0,0,0;
         0,0,0.100, -0.0138,0,0;
         0,0,0.0019, 0.0002, -0.0124,0;
         0,0,0,0.0001, 0.0028, -0.0049]; %1986 Tuesday Lake
[kT,cifT, c2T, cifGT] = contDS(T)

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
RF(9,8) = 6.8257; RF(9,9) = -17.2972; %Panamanian Rainforest
[kRF, cifRF, c2RF, cifGRF] = contDS(RF)

B = boeing_demo('S');
[kB, cifB, c2B, cifGB] = contDS(B,2)

   




































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block Diagonal Matrix Example
rng('default');    
b1 = 10;
b2 = 10;
%b3 = 10;
%b4 = 20;
%b5 = 20;
res = 10000;
%A = [randn(b1)-(b1+(b1/2)*1i)*eye(b1), zeros(b1, b2); zeros(b2, b1), randn(b2)+(b2+(b2/2)*1i)*eye(b2)];
A = [randn(b1)-5*eye(b1), zeros(b1, b2); zeros(b2, b1), randn(b2)+(10+5*1i)*eye(b2)];
%  A = [randn(b1)+(-20)*eye(b1), zeros(b1, b2+b3);%+b4+b5); 
%      zeros(b2, b1), randn(b2), zeros(b2, b3);%+b4+b5); 
%     zeros(b3, b1+b2), randn(b3)+(20)*eye(b3)];%, zeros(b3,b4+b5);
%     zeros(b4, b1+b2+b3), randn(b4)+4*(-20+10i)*eye(b4), zeros(b4, b5);
%     zeros(b5, b1+b2+b3+b4) randn(b5)+4*(20-10i)*eye(b3)];

%%
figure()
[nr, nr_prime] = numerical_range(A, res);
plot(nr), hold on, axis equal
evs = eig(A);
plot(evs, 'kx')
plot(numerical_range(A(1:b1,1:b1), res))
plot(numerical_range(A(b1+1:b1+b2,b1+1:b1+b2), res))
plot(numerical_range(A(b1+b2+1:b1+b2+b3,b1+b2+1:b1+b2+b3), res))

pause
[k, cif, delOm, delOm_prime] = remove_disks(A);
k, cif


%%
%Needs two disks removed to split into two disjoint subsets.
load('tripleBlock0.mat');
%omA = [4+1.5i, 3+4i]; omA2 = [4+2i, 3+4i];

[k, cif, delOmA,~,c1, c2] = remove_disks(A);
c1, c2, k, cif

eigA = eig(A);
plot(eigA, 'rx')

[n, ~] = size(A);
r1 = r_of_A(A, n, omA(1));
r2 = r_of_A(A, n, omA(2));

diskwriter =@(r, om, x) r*exp(1i*x) + om;
d1 = diskwriter(r1, omA(1), linspace(0,2*pi, 100));
d2 = diskwriter(r2, omA(2), linspace(0,2*pi, 100));
%plot disks
plot(d1, 'b--')
hold on
plot(d2, 'b--')

A1 = A(1:10, 1:10); A2 = A(11:end, 11:end);
nr1 = numerical_range(A1,1000);
nr2 = numerical_range(A2,1000);
nr = numerical_range(A, 1000);
plot(nr1, 'r')
plot(nr2, 'r')
plot(nr, 'b--')
plot(delOmA, 'k', 'LineWidth', 1.25)

% Now I calculate analytic formula for when one disk would split the nr
[f1, i1] = max(real(nr1));
[d, i2] = min(abs(nr2-f1));
%min_sig = svds(A, 1, 'smallest');
%xi = nr1(i1) + (nr2(i2) - nr1(i1))/2;
%[fabscissa, i3] = max(abs(nr - xi));
alpha1 = max(imag(nr1))/d;
alpha2 = (max(real(nr2)) - min(real(nr2)))/d;
if(1/(1+2*alpha1) > (alpha1+alpha2)/2)
    test = "one disk should work"
else
    test = "one disk shouldn't work"
end

%% Split a block diagonal matrix into three disjoint sets
load('tripleBlock1.mat')
nr = numerical_range(A,100);
evA = eig(A);
figure()
plot(nr, 'b--')
hold on, axis equal
plot(evA, 'rx')

%Remove the disks at -9.5 and 10 for optimal K value
[k, cif, delOmA, ~, c1, c2] = remove_disks(A);
c1, c2, k, cif

eigA = eig(A);
plot(eigA, 'rx')

[n, ~] = size(A);
r1 = r_of_A(A, n, omA(1));
r2 = r_of_A(A, n, omA(2));

diskwriter =@(r, om, x) r*exp(1i*x) + om;
d1 = diskwriter(r1, omA(1), linspace(0,2*pi, 100));
d2 = diskwriter(r2, omA(2), linspace(0,2*pi, 100));
%plot disks
plot(d1, 'b--')
hold on
plot(d2, 'b--')

A1 = A(1:10, 1:10); A2 = A(11:15, 11:15); A3 = A(16:25, 16:25);
nr1 = numerical_range(A1,1000);
nr2 = numerical_range(A2,1000);
nr3 = numerical_range(A3,1000);
nr  = numerical_range(A, 1000);
plot(nr1, 'r')
plot(nr2, 'r')
plot(nr3, 'r')
plot(nr, 'b--')
plot(delOmA, 'k', 'LineWidth', 1.25)

% between nr1 and nr2
test = "Case 1: Disk containing W(A)"
[f1, i1] = max(real(nr1));
[d, i2] = min(abs(nr2-f1));
d = d/2;
alpha1 = (max(real(nr1)) - min(real(nr1)))/2/d;
alpha2 = (max(real(nr2)) - min(real(nr2)))/2/d;
if(1/(1+2*alpha1) > (alpha1+alpha2)/2)
    test = "one disk should work"
else
    test = "one disk shouldn't work"
end

%between nr2 and nr3
[f1, i1] = max(real(nr2));
[d, i2] = min(abs(nr3-f1));
d = d/2;
alpha1 = (max(real(nr2)) - min(real(nr2)))/2/d;
alpha2 = (max(real(nr3)) - min(real(nr3)))/2/d;
if(1/(1+2*alpha1) > (alpha1+alpha2)/2)
    test = "one disk should work"
else
    test = "one disk shouldn't work"
end

% between nr1 and nr2 min "spanning" disk instead
test = "Case 2: min Disk tangent to W(A) at least twice."
[f1, i1] = max(real(nr1));
[d, i2] = min(abs(nr2-f1));
d = d/2;
alpha1 = (max(imag(nr1)))/2/d;
alpha2 = (max(imag(nr2)))/2/d;
if(1/(1+2*alpha1) > (alpha1+alpha2)/2)
    test = "one disk should work"
else
    test = "one disk shouldn't work"
end

%between nr2 and nr3
[f1, i1] = max(real(nr2));
[d, i2] = min(abs(nr3-f1));
d = d/2;
alpha1 = (max(imag(nr2)))/2/d;
alpha2 = (max(imag(nr3)))/2/d;
if(1/(1+2*alpha1) > (alpha1+alpha2)/2)
    test = "one disk should work"
else
    test = "one disk shouldn't work"
end

%Maybe the min "spanning" disk is a batter way of thinking of it rather
%than a disk containing W(A)!

%% Single Disk Example

rng('default')
B11 = randn(4);
B22 = randn(4) + 8*eye(4);
B = [B11, zeros(4); zeros(4), B22];

eigB = eig(B);

[n, ~] = size(B);

nr1 = numerical_range(B11,1000);
nr2 = numerical_range(B22,1000);
[nr, nrprime]  = numerical_range(B, 1000);
[nr, nrprime] = nr_gap_fill(nr, nrprime);

zeta = 3.5;
[~, r1] = r_of_A(B, n, zeta);
diskwriter =@(r, om, x) r*exp(1i*x) + om;
d1 = diskwriter(r1, zeta, linspace(0,2*pi, 100));

[delOmB, delomvec, xs, r1orr2] = define_del_Omega({nr}, {zeros(1, length(nr))}, B, zeta, 1000, r1);
[delOmB, delomvec] = cellmat2plot(delOmB, delomvec,1);
[k, cif, delOm_prime, c1, c2] = calc_k_removed_disks(B, zeta, nr, nrprime, delOmB, delomvec, xs, r1orr2);
k, cif, c1, c2

%plot disks
plot(d1, 'b--')
hold on
%plot eigs
plot(eigB, 'rx')
%plot numerical ranges
plot(nr1, 'r')
plot(nr2, 'r')
plot(nr, 'b--')
plot(delOmB, 'k', 'LineWidth', 1.25)
axis equal






%% Gracar 100 Example
n = 100;
A = gallery('grcar',n);

eigA = eig(A);

[n, ~] = size(A);

[nr, nrprime]  = numerical_range(A, 1000);
[nr, nrprime] = nr_gap_fill(nr, nrprime);

zeta = 0;
[~, r1] = r_of_A(A, n, zeta);
diskwriter =@(r, zeta, x) r*exp(1i*x) + zeta;
d1 = diskwriter(r1, zeta, linspace(0,2*pi, 100));

[delOmA, delomvec, xs, r1orr2] = define_del_Omega({nr}, {zeros(1, length(nr))}, A, zeta, 10000, r1);
[delOmA, delomvec] = cellmat2plot(delOmA, delomvec,1);
[k, cif, delOm_prime] = calc_k_removed_disks(A, zeta, nr, nrprime, delOmA, delomvec, xs, r1orr2);
k, cif

%plot disks
plot(d1, 'b--')
hold on
%plot eigs
plot(eigA, 'rx')
%plot numerical ranges
plot(nr, 'b--')
plot(delOmA, 'k', 'LineWidth', 1.25)
axis equal

%% Gracar 32 Example
n = 32;
%For Grchar(32) use epsilon = 10^-3
%For Grchar(24) use epsilon = 10^-2.2
eps = -3; %log 10 of the epsilon value.
A = gallery('grcar',n);

eigA = eig(A);

% opts.levels = [-3]; opts.ax = [-0.5 2 -2.6 2.6]; opts.npts = 2000;
% eigtool(A,opts)
% pause
% GamA = pe_contour(xGrchar, yGrchar, ZGrchar, 10.^[eps, eps], 1);
% delOmA = GamA(2:2:end);
% delOmA = cellmat2plot(delOmA,1);
% delOmA = delOmega_flipper(delOmA,1);
% delOmAP = delOmprime2(delOmA);
% save('psGrcharBIG.mat', 'delOmA', 'delOmAP', 'xGrchar','yGrchar','ZGrchar')

load('psGrcharBIG.mat')

%calculate the two different values of K for the set of epsilon=10^-3 pseudospectra 
[c1A, c1_index] = calc_c1(delOmA, delOmAP)
breaks = find(isnan(delOmA)); breaks = [0,breaks, length(delOmA)+1];
c2A = 0; cifA = 0; L = 0;
for jj = 1:length(breaks)-1
    c2A = c2A + calc_c2_curve(A, delOmA(breaks(jj)+1:breaks(jj+1)-1), delOmAP(breaks(jj)+1:breaks(jj+1)-1));
    cifA = cifA + resolvent_norm_integral(A, delOmA(breaks(jj)+1:breaks(jj+1)-1));
    L = L + sum(abs([delOmA(breaks(jj)+2:breaks(jj+1)-1), delOmA(breaks(jj)+1)]-delOmA(breaks(jj)+1:breaks(jj+1)-1)));
end
c2A
kA = c2A + sqrt(c1A + c2A^2)
cifA
cifAexact = L/(2*pi)*10^-eps










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pseudospectral Example

%The goal is to use the epsilon pseudospectra to define Omega for a matrix.
% Then the cif has a closed form analytic expression, but Omega is no
% longer convex.

A = [0 1 2; -0.01 0 3; 0 0 0]; %eps = 10^-4
% opts.levels = [-4]; opts.ax = [-0.03 0.03 -0.12 0.12]; opts.npts = 1000;
% eigtool(A, opts)
% GamA = pe_contour(xA, yA, ZA, 10.^[-4, -4], 1);
% delOmA = GamA(2:2:end);
% delOmA = cellmat2plot(delOmA,1);
% delOmA = delOmega_flipper(delOmA,1);
B = [-5 4 4; 0 -2-2i 4; 0 0 -0.3+1i]; %eps = 10^-0.35
% opts.levels = [-0.35]; opts.ax = [-5.7, 0.6 -2.9 1.9]; opts.npts = 1000;
% eigtool(B, opts)
% GamB = pe_contour(xB, yB, ZB, 10.^[-0.35, -0.35], 1);
% delOmB = GamB(2:2:end);
% delOmB = cellmat2plot(delOmB,1);
% delOmB = delOmega_flipper(delOmB, 1);
% 
% delOmAP = delOmprime2(delOmA);
% delOmBP = delOmprime2(delOmB);
% save('psOmegas.mat', 'delOmB', 'delOmA', 'delOmAP', 'delOmBP')

load('psOmegas.mat')

% For matrix A and Omega equal to the 10^-3.9 pseudospectral set find
% K
c1A = calc_c1(delOmA, delOmAP);
breaks = find(isnan(delOmA));
c2A = calc_c2_curve(A, delOmA(1:breaks(1)-1), delOmAP(1:breaks(1)-1))+...
    calc_c2_curve(A, delOmA(breaks(1)+1:breaks(2)-1), delOmAP(breaks(1)+1:breaks(2)-1))+...
    calc_c2_curve(A, delOmA(breaks(2)+1:end), delOmAP(breaks(2)+1:end));
kA = c2A + sqrt(c1A + c2A^2)
% the integral of the resolvent norm
cifA = resolvent_norm_integral(A, delOmA(1:breaks(1)-1))+...
    resolvent_norm_integral(A, delOmA(breaks(1)+1:breaks(2)-1))+...
    resolvent_norm_integral(A, delOmA(breaks(2)+1:end))
%analytic formula for the integral of the resolvent norm 
L = cat(2, [delOmA(2:breaks(1)-1), delOmA(1)] - delOmA(1:breaks(1)-1),...
    [delOmA(breaks(1)+2:breaks(2)-1), delOmA(breaks(1)+1)]-delOmA(breaks(1)+1:breaks(2)-1), ...
    [delOmA(breaks(2)+2:end), delOmA(breaks(2)+1)]-delOmA(breaks(2)+1:end));
L = sum(abs(L));
cifAexact = L/(2*pi)*10^4

% For matrix B and Omega equal to the 10^-0.35 pseudospectral set find
% K
c1B = calc_c1(delOmB, delOmBP);
breaks = find(isnan(delOmB));
c2B = calc_c2_curve(B, delOmB(1:breaks(1)-1), delOmBP(1:breaks(1)-1))+...
    calc_c2_curve(B, delOmB(breaks(1)+1:end), delOmBP(breaks(1)+1:end));
kB = c2B + sqrt(c1B + c2B^2)
% the integral of the resolvent norm
cifB = resolvent_norm_integral(B, delOmB(1:breaks(1)-1))+...
    resolvent_norm_integral(B, delOmB(breaks(1)+1:end))
%analytic formula for the integral of the resolvent norm 
L = cat(2, [delOmB(2:breaks(1)-1), delOmB(1)] - delOmB(1:breaks(1)-1),...
    [delOmB(breaks(1)+2:end), delOmB(breaks(1)+1)]-delOmB(breaks(1)+1:end));
L = sum(abs(L));
cifBexact = L/(2*pi)*10^0.35


%% Tuesday Lake with different types of Omega sets

T = [-0.9503,0,0.0690,0.0002, 0.0027,0.0034;
         0.95, -0.18,0,0,0,0;
         0, 0.15, -0.2569,0,0,0;
         0,0,0.100, -0.0138,0,0;
         0,0,0.0019, 0.0002, -0.0124,0;
         0,0,0,0.0001, 0.0028, -0.0049]; %1986 Tuesday Lake
[kT,cifT] = contDS(T)

% opts.levels = [-1.1, -1.3] ; opts.ax = [-1.06 0.08 -0.14 0.14] ; opts.npts = 1000;
% eigtool(T, opts)
% export the pseudospectral data from eigtool
load("TuesdayPseudospectra.mat")
GamT1 = pe_contour(xT, yT, ZT, 10.^[-1.1, -1.1], 1);
GamT2 = pe_contour(xT2,yT2,ZT2, 10.^[-2.35,-2.35],1);
GamT3 = pe_contour(xT, yT, ZT, 10.^[-1.3, -1.3], 1);
delOmT1 = GamT1(2:2:end);
delOmT1 = cellmat2plot(delOmT1, 1);
delOmT1 = delOm_flipper(delOmT1,1);
delOmT1_prime = calc_delOm_prime(delOmT1);
delOmT2 = GamT2(2:2:end);
delOmT2 = cellmat2plot(delOmT2, 1);
delOmT2 = delOm_flipper(delOmT2,1);
delOmT2_prime = calc_delOm_prime(delOmT2);
delOmT3 = GamT3(2:2:end);
delOmT3 = cellmat2plot(delOmT3, 1);
delOmT3 = delOm_flipper(delOmT3,1);
delOmT3_prime = calc_delOm_prime(delOmT3);

%epsilon = -1.1
c1T1 = calc_c1(delOmT1, delOmT1_prime);
breaks = find(isnan(delOmT1));
%technically this set does extend past nrT, so this is a overestimation of c2
c2T1 = calc_c2_curve(T, delOmT1(1:breaks(1)-1), delOmT1_prime(1:breaks(1)-1))+...
    calc_c2_curve(T, delOmT1(breaks(1)+1:end), delOmT1_prime(breaks(1)+1:end));
kT1 = c2T1 + sqrt(c1T1 + c2T1^2)
% the integral of the resolvent norm
cifT1 = resolvent_norm_integral(T, delOmT1(1:breaks(1)-1))+...
    resolvent_norm_integral(T, delOmT1(breaks(1)+1:end))
%analytic formula for the integral of the resolvent norm 
L = cat(2, [delOmT1(2:breaks(1)-1), delOmT1(1)] - delOmT1(1:breaks(1)-1),...
    [delOmT1(breaks(1)+2:end), delOmT1(breaks(1)+1)]-delOmT1(breaks(1)+1:end));
L = sum(abs(L));
cifT1exact = L/(2*pi)*10^1.1


%epsilon = -2.35
%For a set that remains in the left-half plane
c1T2 = calc_c1(delOmT2, delOmT2_prime);
breaks = find(isnan(delOmT2));
%estimate c2
c2T2 = calc_c2_curve(T, delOmT2(1:breaks(1)-1), delOmT2_prime(1:breaks(1)-1))+...
    calc_c2_curve(T, delOmT2(breaks(1)+1:breaks(2)-1), delOmT2_prime(breaks(1)+1:breaks(2)-1))+...
    calc_c2_curve(T, delOmT2(breaks(2)+1:breaks(3)-1), delOmT2_prime(breaks(2)+1:breaks(3)-1))+...
    calc_c2_curve(T, delOmT2(breaks(3)+1:end), delOmT2_prime(breaks(3)+1:end));
kT2 = c2T2 + sqrt(c1T2 + c2T2^2)
% the integral of the resolvent norm
cifT2 = resolvent_norm_integral(T, delOmT1(1:breaks(1)-1))+...
    resolvent_norm_integral(T, delOmT2(breaks(1)+1:breaks(2)-1))+...
    resolvent_norm_integral(T, delOmT2(breaks(2)+1:breaks(3)-1))+...
    resolvent_norm_integral(T, delOmT2(breaks(3)+1:end))
%analytic formula for the integral of the resolvent norm 
% L = cat(2, [delOmT1(2:breaks(1)-1), delOmT1(1)] - delOmT1(1:breaks(1)-1),...
%     [delOmT1(breaks(1)+2:end), delOmT1(breaks(1)+1)]-delOmT1(breaks(1)+1:end));
% L = sum(abs(L));
% cifT1exact = L/(2*pi)*10^1.1

%epsilon = -1.3
c1T3 = calc_c1(delOmT3, delOmT3_prime);
breaks = find(isnan(delOmT3));
c2T3 = calc_c2_curve(T, delOmT3(1:breaks(1)-1), delOmT3_prime(1:breaks(1)-1))+...
    calc_c2_curve(T, delOmT3(breaks(1)+1:breaks(2)-1), delOmT3_prime(breaks(1)+1:breaks(2)-1))+...
    calc_c2_curve(T, delOmT3(breaks(2)+1:end), delOmT3_prime(breaks(2)+1:end));
kT3 = c2T3 + sqrt(c1T3 + c2T3^2)
% the integral of the resolvent norm
cifT3 = resolvent_norm_integral(T, delOmT3(1:breaks(1)-1))+...
    resolvent_norm_integral(T, delOmT3(breaks(1)+1:breaks(2)-1))+...
    resolvent_norm_integral(T, delOmT3(breaks(2)+1:end))
%analytic formula for the integral of the resolvent norm 
L = cat(2, [delOmT3(2:breaks(1)-1), delOmT3(1)] - delOmT3(1:breaks(1)-1),...
    [delOmT3(breaks(1)+2:breaks(2)-1), delOmT3(breaks(1)+1)]-delOmT3(breaks(1)+1:breaks(2)-1), ...
    [delOmT3(breaks(2)+2:end), delOmT3(breaks(2)+1)]-delOmT3(breaks(2)+1:end));
L = sum(abs(L));
cifT3exact = L/(2*pi)*10^1.3




%%
load("transient9_pseudo.mat");

GamDem = pe_contour(xDem, yDem, ZDem, 10.^[-1.55, -1.55], 0);
delOmDem = cat(2, GamDem{1,2}, NaN, GamDem{1,4});
breaks = find(isnan(delOmDem));

L = cat(2, [delOmDem(2:breaks-1), delOmDem(1)] - delOmDem(1:breaks-1),...
    [delOmDem(breaks+2:end), delOmDem(breaks+1)]-delOmDem(breaks+1:end));
L = sum(abs(L));
cifDemexact = L/(2*pi)*10^1.55

































%% Advection with Dirichlet boundary conditions only on one side
% note with periodic boundary conditions, this problem results in a normal
% matrix. With only one boundary given, then it results in a non-normal
% matrix

load("advectionEx.mat");
%% To rerun the code instead of loading the experiment

%define the matrix; n (integer) is the dimension aka number of steps
Aad = @(n) diag(zeros(n,1)) + diag(ones(n-1,1),1) + diag(-1*ones(n-1,1),-1);

dim = 7;
k = zeros(1,dim); cif = zeros(1,dim); 
for n = 3:10
    A = Aad(n);
    A(end, end-2:end) = [1 -4 3];
    A = A - 2*eye(n);

    [k(n-2), cif(n-2)] = contDS(A);
end
k
cif


dim = 3;
kbig = zeros(1,dim); cifbig = zeros(1,dim); 
for n = 50:50:150
    A = Aad(n);
    A(end, end-2:end) = [1 -4 3];
    A = A - 2*eye(n);

    [kbig(n/50), cifbig(n/50)] = contDS(A);
end
kbig
cifbig

save("advectionEx.mat", "n", "k", "cif", "nbig", "kbig", "cifbig")


%%
P = [.5 .5 0; .125 .75 .125; 0 .5 .5];
% opts.levels = [-1.2]; opts.ax = [.15 1.1 -0.1 0.1]; opts.npts = 1000;
% eigtool(P, opts)
% GamP = pe_contour(xP, yP, ZP, 10.^[-1.2, -1.2], 1);
% delOmP = GamP(2:2:end);
% delOmP = cellmat2plot(delOmP,1);
% delOmP = delOmega_flipper(delOmP,1);
% delOmP_prime = delOmprime2(delOmP);

% For matrix P and Omega equal to the 10^-1.2 pseudospectral set find
% K
c1P = calc_c1(delOmP, delOmP_prime);
breaks = find(isnan(delOmP));
c2P = calc_c2_curve(P, delOmP(1:breaks(1)-1), delOmP_prime(1:breaks(1)-1))+...
    calc_c2_curve(P, delOmP(breaks(1)+1:breaks(2)-1), delOmP_prime(breaks(1)+1:breaks(2)-1))+...
    calc_c2_curve(P, delOmP(breaks(2)+1:end), delOmP_prime(breaks(2)+1:end));
kP = c2P + sqrt(c1P + c2P^2)
% the integral of the resolvent norm
cifP = cauchyIntFormula(P, delOmP(1:breaks(1)-1))+...
    cauchyIntFormula(P, delOmP(breaks(1)+1:breaks(2)-1))+...
    cauchyIntFormula(P, delOmP(breaks(2)+1:end))
%analytic formula for the integral of the resolvent norm 
L = cat(2, [delOmP(2:breaks(1)-1), delOmP(1)] - delOmP(1:breaks(1)-1),...
    [delOmP(breaks(1)+2:breaks(2)-1), delOmP(breaks(1)+1)]-delOmP(breaks(1)+1:breaks(2)-1), ...
    [delOmP(breaks(2)+2:end), delOmP(breaks(2)+1)]-delOmP(breaks(2)+1:end));
L = sum(abs(L));
cifPexact = L/(2*pi)*10^1.2