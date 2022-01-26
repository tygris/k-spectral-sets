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
Gam_tdl = cell(1,11); td_lbl = cell(1,11); td_upl = cell(1,11);
dims = 10:4:50;
for jj = 1:6
    Gam_td = pe_contour(cell2mat(x_td(jj)),cell2mat(y_td(jj)),cell2mat(Z_td(jj)),cell2mat(powertd(jj)), 0);
    td_lb = pseudo_lb(Gam_td, 'c');
    td_up = exp(1)*dims(jj)*max(td_lb(2, :));
    Gam_tdl(jj) = {Gam_td};
    td_lbl(jj) = {td_lb};
    td_upl(jj) = {td_up};
end

% Now find our upper bound on the matrix envelope
td_kl = cell(1,11); td_cifl = cell(1,11); td_c2l = cell(1,11);
for jj = 1:6
    A = transient_demo(4*jj+6);
    %calculate spectral set
    [k, cif, c2] = contDS(A,4); 
    td_kl(jj) = {k};
    td_cifl(jj) = {cif};
    td_c2l(jj) = {c2};
end
    
td_upl
td_kl
td_cifl
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Using the contDS Function pass these matrices

A = [.5 -1; 1, -1.5];
[kA, cifA] = contDS(A)

C = [0 1 2; -0.01 0 3; 0 0 0];
[kC, cifC] = contDS(C)

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

%set-up for removing a disk
[nr, nr_prime] = nrGapFill(nr, nr_prime);
delOm = {nr};
delom = {zeros(1,length(nr))};

mat = A;
res_num = 2000;

[numRange, nr_prime] = numerical_range(mat,res_num);
[numRange, nr_prime] = nrGapFill(numRange, nr_prime);
figure()
plot(numRange), daspect([1,1,1])


%% 2. The user is asked where they would like the center of removed disks to be

Y = 'Y'; N = 'N';
om = [];
r_over_pi = [];
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
        [del_Om, del_om, xs_new, roverpi_new] = define_del_Omega(del_Om, del_om, mat, om_new, res_num);
        om = cat(2, om, om_new);
        r_over_pi = cat(2, r_over_pi, roverpi_new);
        xs = cat(2, xs, xs_new);
    end
    more = input('Would you like to remove more disks? (Y/N)\n');
    if more ~= 'Y'
        if 'N' == input('Are you sure? (Y/N)\n')
            more = input('Would you like to remove more disks? (Y/N)\n');
        end
    end
end

figure()
plot(numRange, '--k'), hold on, axis equal
plot(cellmat2plot(del_Om,1),'b')
