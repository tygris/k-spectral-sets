%File of different examples of how K-Spectral Sets can be used

%Natalie Wellen
%1/25/21

%% Example 1
% Using the transient_demo() from Eigtool: https://www.cs.ox.ac.uk/pseudospectra/eigtool/

A =  transient_demo(10);

[k, cif] = contDS(A,2)

[k, cif] = discDS(A, .75)

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
    figure()
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
dims = 10:4:50; dimlength = length(dims);
Gam_tdl = cell(1,dimlength); td_lbl = cell(1,dimlength); td_upl = zeros(1,dimlength);
for jj = 1:dimlength
    Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :));
    Gam_tdl(jj) = {Gam_td};
    td_lbl(jj) = {td_lb};
    td_upl(jj) = td_up;
end

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
semilogy(dims, [td_ups td_upl], 'DisplayName', 'Pseudospectral', opts{:})
hold on
plot(dims, [td_ks td_kl], 'DisplayName', 'K', opts{:})
plot(dims, [td_cifs td_cifl], 'DisplayName', 'Resolvent Norm', opts{:})
title('Transient Matrix Growth Bounds')
xlabel('dim(A)')
ylabel('Upper Bound of ||f(A)||')
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Using the contDS Function pass these matrices

A = [.5 -1; 1, -1.5];
[kA, cifA] = contDS(A)

C = [0 1 2; -0.01 0 3; 0 0 0];
[kC, cifC] = discDS(C)
%C is a bad example for what we do because one of its eigenvalues = 0+0i

T = [-0.9503,0,0.0690,0.0002, 0.0027,0.0034;
         0.95, -0.18,0,0,0,0;
         0, 0.15, -0.2569,0,0,0;
         0,0,0.100, -0.0138,0,0;
         0,0,0.0019, 0.0002, -0.0124,0;
         0,0,0,0.0001, 0.0028, -0.0049]; %1986 Tuesday Lake
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
RF(9,8) = 6.8257; RF(9,9) = -17.2972; %Panamanian Rainforest
[kRF, cifRF] = contDS(RF)

B = boeing_demo('S');
[kB, cifB] = contDS(B,2)

   




































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block Diagonal Matrix Example
    
b1 = 10;
b2 = 10;
%b3 = 2;
res = 10000;
A = [randn(b1)-5*eye(b1), zeros(b1, b2); zeros(b2, b1), randn(b2)+(10+5i)*eye(b2)];
%A = [randn(b1), zeros(b1, b2+b3); zeros(b2, b1), randn(b2)+(10+5i)*eye(b2), zeros(b2, b3); zeros(b3, b1+b2), rand(b3)-(10-10*1i)*eye(b3)];


figure()
[nr, nr_prime] = numerical_range(A, res);
plot(nr), hold on
evs = eig(A);
plot(evs, 'kx')
plot(numerical_range(A(1:b1,1:b1), res))
plot(numerical_range(A(b1+1:b1+b2,b1+1:b1+b2), res))
%plot(numerical_range(A(b1+b2+1:b1+b2+b3,b1+b2+1:b1+b2+b3), res))

%% Remove disks from A
% [nr, nr_prime] = nrGapFill(nr, nr_prime);
% delOm = {nr};
% delom = {zeros(1,length(nr))};

mat = A;
res_num = 2000;

[numRange, nr_prime] = numerical_range(mat,res_num);
[numRange, nr_prime] = nrGapFill(numRange, nr_prime);
figure()
plot(numRange), daspect([1,1,1])


% 2. The user is asked where they would like the center of removed disks to be

Y = 'Y'; N = 'N';
om = [];
r1orr2 = [];
xs = [];
more = 'Y';
del_Om = {numRange};
del_om = {zeros(1,length(numRange))};
moveon = 0;
while more == 'Y'
    om_new = input("Where would you like to remove a disk(s)?\n");
    close
    del_Om_vec = cellmat2plot(del_Om,1);
    figure()
    plot(real(del_Om_vec), imag(del_Om_vec))
    daspect([1,1,1]);
    hold on
    plot(real(om_new), imag(om_new), 'mo');
    moveon = input('Is this where you would like to remove the disk? (Y/N)\n');
    if moveon == 'Y'
        close
        [del_Om, del_om, xs_new, r1orr2_new] = define_del_Omega(del_Om, del_om, mat, om_new, res_num);
        om = cat(2, om, om_new);
        r1orr2 = cat(2, r1orr2, r1orr2_new);
        xs = cat(2, xs, xs_new);
    end
    more = input('Would you like to remove more disks? (Y/N)\n');
    if more ~= 'Y'
        if 'N' == input('Are you sure? (Y/N)\n')
            more = input('Would you like to remove more disks? (Y/N)\n');
        end
    end
end

%define delOm
[delOmvec, delomvec] = cellmat2plot(del_Om, del_om);
figure()
plot(numRange, '--k'), hold on, axis equal
plot(delOmvec,'b')

%calculate K and the integral of the resolvent norm for A
%define delOmPrime
indBreaks = find(ismember(numRange,xs));
delOm_prime = zeros(1, length(delOmvec));
    %first along points of delOm that coincide with the numerical range
primer = [];
if ~delomvec(1) %false implies 1:indBreaks(1) is nr_prime
    indBreaks = cat(2, 1, indBreaks, length(numRange));
end
for jj = 1:2:length(indBreaks)-1
    primer = cat(2, primer, nr_prime(indBreaks(jj):indBreaks(jj+1)));
end
delOm_prime(delomvec == 0) = primer;
    %then with points along boundaries of the removed disks
for jj = 1:max(delomvec)
    om_now = om(jj);
    kk = find(delomvec == jj);
    radius = abs(delOmvec(kk(1)) - om_now);
    delOm_prime(kk) = -1i*((delOmvec(kk) - om_now)/radius);
end
% c1 is numerically estimated
[c1] = calc_c1(delOmvec, delOm_prime);
%We use Thoerem 2 to calculate a bound on c2
c2 = 1+sum(r1orr2);
k = c2+sqrt(c1+c2^2)
indBreaks = find(isnan(delOmvec));
indBreaks = cat(2, indBreaks, length(delOmvec)+1);
cif = cauchyIntFormula(mat, delOmvec(1:indBreaks(1)-1));
for jj = 2:length(indBreaks)
    cif = cif + cauchyIntFormula(mat, delOmvec(indBreaks(jj-1)+1:indBreaks(jj)-1));
end
cif


%%
%Needs two disks removed.
A =   [-4.4623e+00 + 0.0000e+00i	-1.3499e+00 + 0.0000e+00i	6.7150e-01 + 0.0000e+00i	8.8840e-01 + 0.0000e+00i	-1.0224e-01 + 0.0000e+00i	-8.6365e-01 + 0.0000e+00i	-1.0891e+00 + 0.0000e+00i	-6.1560e-01 + 0.0000e+00i	1.4193e+00 + 0.0000e+00i	-1.1480e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
1.8339e+00 + 0.0000e+00i	-1.9651e+00 + 0.0000e+00i	-1.2075e+00 + 0.0000e+00i	-1.1471e+00 + 0.0000e+00i	-2.4145e-01 + 0.0000e+00i	7.7359e-02 + 0.0000e+00i	3.2557e-02 + 0.0000e+00i	7.4808e-01 + 0.0000e+00i	2.9158e-01 + 0.0000e+00i	1.0487e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-2.2588e+00 + 0.0000e+00i	7.2540e-01 + 0.0000e+00i	-4.2828e+00 + 0.0000e+00i	-1.0689e+00 + 0.0000e+00i	3.1921e-01 + 0.0000e+00i	-1.2141e+00 + 0.0000e+00i	5.5253e-01 + 0.0000e+00i	-1.9242e-01 + 0.0000e+00i	1.9781e-01 + 0.0000e+00i	7.2225e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
8.6217e-01 + 0.0000e+00i	-6.3055e-02 + 0.0000e+00i	1.6302e+00 + 0.0000e+00i	-5.8095e+00 + 0.0000e+00i	3.1286e-01 + 0.0000e+00i	-1.1135e+00 + 0.0000e+00i	1.1006e+00 + 0.0000e+00i	8.8861e-01 + 0.0000e+00i	1.5877e+00 + 0.0000e+00i	2.5855e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.1877e-01 + 0.0000e+00i	7.1474e-01 + 0.0000e+00i	4.8889e-01 + 0.0000e+00i	-2.9443e+00 + 0.0000e+00i	-5.8649e+00 + 0.0000e+00i	-6.8493e-03 + 0.0000e+00i	1.5442e+00 + 0.0000e+00i	-7.6485e-01 + 0.0000e+00i	-8.0447e-01 + 0.0000e+00i	-6.6689e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-1.3077e+00 + 0.0000e+00i	-2.0497e-01 + 0.0000e+00i	1.0347e+00 + 0.0000e+00i	1.4384e+00 + 0.0000e+00i	-3.0051e-02 + 0.0000e+00i	-3.4674e+00 + 0.0000e+00i	8.5931e-02 + 0.0000e+00i	-1.4023e+00 + 0.0000e+00i	6.9662e-01 + 0.0000e+00i	1.8733e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-4.3359e-01 + 0.0000e+00i	-1.2414e-01 + 0.0000e+00i	7.2689e-01 + 0.0000e+00i	3.2519e-01 + 0.0000e+00i	-1.6488e-01 + 0.0000e+00i	-7.6967e-01 + 0.0000e+00i	-6.4916e+00 + 0.0000e+00i	-1.4224e+00 + 0.0000e+00i	8.3509e-01 + 0.0000e+00i	-8.2494e-02 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.4262e-01 + 0.0000e+00i	1.4897e+00 + 0.0000e+00i	-3.0344e-01 + 0.0000e+00i	-7.5493e-01 + 0.0000e+00i	6.2771e-01 + 0.0000e+00i	3.7138e-01 + 0.0000e+00i	-7.4230e-01 + 0.0000e+00i	-4.5118e+00 + 0.0000e+00i	-2.4372e-01 + 0.0000e+00i	-1.9330e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.5784e+00 + 0.0000e+00i	1.4090e+00 + 0.0000e+00i	2.9387e-01 + 0.0000e+00i	1.3703e+00 + 0.0000e+00i	1.0933e+00 + 0.0000e+00i	-2.2558e-01 + 0.0000e+00i	-1.0616e+00 + 0.0000e+00i	-1.7738e-01 + 0.0000e+00i	-4.7843e+00 + 0.0000e+00i	-4.3897e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
2.7694e+00 + 0.0000e+00i	1.4172e+00 + 0.0000e+00i	-7.8728e-01 + 0.0000e+00i	-1.7115e+00 + 0.0000e+00i	1.1093e+00 + 0.0000e+00i	1.1174e+00 + 0.0000e+00i	2.3505e+00 + 0.0000e+00i	-1.9605e-01 + 0.0000e+00i	-1.1658e+00 + 0.0000e+00i	-6.7947e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.0840e+01 + 5.0000e+00i	-2.1384e+00 + 0.0000e+00i	2.9080e+00 + 0.0000e+00i	-3.5385e-01 + 0.0000e+00i	2.2890e-02 + 0.0000e+00i	5.2006e-01 + 0.0000e+00i	-2.9375e-01 + 0.0000e+00i	-1.3320e+00 + 0.0000e+00i	-1.3617e+00 + 0.0000e+00i	-1.9522e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-8.8803e-01 + 0.0000e+00i	9.1604e+00 + 5.0000e+00i	8.2522e-01 + 0.0000e+00i	-8.2359e-01 + 0.0000e+00i	-2.6200e-01 + 0.0000e+00i	-2.0028e-02 + 0.0000e+00i	-8.4793e-01 + 0.0000e+00i	-2.3299e+00 + 0.0000e+00i	4.5503e-01 + 0.0000e+00i	-2.1761e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.0009e-01 + 0.0000e+00i	1.3546e+00 + 0.0000e+00i	1.1379e+01 + 5.0000e+00i	-1.5771e+00 + 0.0000e+00i	-1.7502e+00 + 0.0000e+00i	-3.4771e-02 + 0.0000e+00i	-1.1201e+00 + 0.0000e+00i	-1.4491e+00 + 0.0000e+00i	-8.4871e-01 + 0.0000e+00i	-3.0311e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-5.4453e-01 + 0.0000e+00i	-1.0722e+00 + 0.0000e+00i	-1.0582e+00 + 0.0000e+00i	1.0508e+01 + 5.0000e+00i	-2.8565e-01 + 0.0000e+00i	-7.9816e-01 + 0.0000e+00i	2.5260e+00 + 0.0000e+00i	3.3351e-01 + 0.0000e+00i	-3.3489e-01 + 0.0000e+00i	2.3046e-02 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	3.0352e-01 + 0.0000e+00i	9.6095e-01 + 0.0000e+00i	-4.6862e-01 + 0.0000e+00i	2.8198e-01 + 0.0000e+00i	9.1686e+00 + 5.0000e+00i	1.0187e+00 + 0.0000e+00i	1.6555e+00 + 0.0000e+00i	3.9135e-01 + 0.0000e+00i	5.5278e-01 + 0.0000e+00i	5.1290e-02 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-6.0033e-01 + 0.0000e+00i	1.2405e-01 + 0.0000e+00i	-2.7247e-01 + 0.0000e+00i	3.3480e-02 + 0.0000e+00i	-9.7921e-01 + 0.0000e+00i	9.8668e+00 + 5.0000e+00i	3.0754e-01 + 0.0000e+00i	4.5168e-01 + 0.0000e+00i	1.0391e+00 + 0.0000e+00i	8.2606e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	4.8997e-01 + 0.0000e+00i	1.4367e+00 + 0.0000e+00i	1.0984e+00 + 0.0000e+00i	-1.3337e+00 + 0.0000e+00i	-1.1564e+00 + 0.0000e+00i	-7.1453e-01 + 0.0000e+00i	8.7429e+00 + 5.0000e+00i	-1.3028e-01 + 0.0000e+00i	-1.1176e+00 + 0.0000e+00i	1.5270e+00 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	7.3936e-01 + 0.0000e+00i	-1.9609e+00 + 0.0000e+00i	-2.7787e-01 + 0.0000e+00i	1.1275e+00 + 0.0000e+00i	-5.3356e-01 + 0.0000e+00i	1.3514e+00 + 0.0000e+00i	-8.6547e-01 + 0.0000e+00i	1.0184e+01 + 5.0000e+00i	1.2607e+00 + 0.0000e+00i	4.6691e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.7119e+00 + 0.0000e+00i	-1.9770e-01 + 0.0000e+00i	7.0154e-01 + 0.0000e+00i	3.5018e-01 + 0.0000e+00i	-2.0026e+00 + 0.0000e+00i	-2.2477e-01 + 0.0000e+00i	-1.7653e-01 + 0.0000e+00i	-4.7615e-01 + 0.0000e+00i	1.0660e+01 + 5.0000e+00i	-2.0971e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-1.9412e-01 + 0.0000e+00i	-1.2078e+00 + 0.0000e+00i	-2.0518e+00 + 0.0000e+00i	-2.9907e-01 + 0.0000e+00i	9.6423e-01 + 0.0000e+00i	-5.8903e-01 + 0.0000e+00i	7.9142e-01 + 0.0000e+00i	8.6202e-01 + 0.0000e+00i	-6.7866e-02 + 0.0000e+00i	1.0625e+01 + 5.0000e+00i];
omA = [4+2i, 3+4i];
