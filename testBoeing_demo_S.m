%script to run for finding the resolvent norm integral and value of K for
%the continuous time DS of the boeing_demo('S')

%1/19/22

A = boeing_demo('S');
res = 20001;
%calculate the numerical range and it's derivative
[nr, nr_prime] = numerical_range(A,res);
%find Gamma1 and it's derivative
[y1, y2] = nr_cut_off(nr, 0); %find the values on the border of the numerical range Gamma1 runs between
Gam1 = flip(chebpts(res, [y2, y1]).');
Gam1 = [y1, Gam1((res-1)/2+2:res-1)+y1, 0, Gam1(2:(res-1)/2)+y2, y2];
Gam1_prime = 1i*ones(1,res);
ds = abs(Gam1(2:res)-Gam1(1:res-1));
figure()
plot(ds)
%define delOmega and it's derivative
ind1 = real(nr)<0;
delOm = cat(2, Gam1, nr(ind1));
%delOm_prime = cat(2, Gam1_prime, nr_prime(ind1));

%calculate c2
%gamma(s) and resolvent norm at each point of Gamma_1
m = length(A);
gammas = zeros(1,res);
rnorms = zeros(1,res);
for jj = 1:res
    R = Gam1_prime(jj)*inv(Gam1(jj)*eye(m) - A); %the matrix s'(sI-A)^-1
    gammas(jj) = abs(min(eig(1/(2*pi*1i)*(R-R')))); 
    rnorms(jj) = norm(R, 2); %/2*pi
end


%Use the trapezoidal rule to estimate the integral of gamma(s)
    %note we assume points are equidistant
%ds = abs(Gam1(2:n)-Gam1(1:n-1));
midpoints = gammas(2:res) + gammas(1:res-1); %/2
midpoints2 = rnorms(2:res) + rnorms(1:res-1); %/2
integral = (sum(midpoints.*ds) +...
    ds(1)*gammas(1)+ds(res-1)*gammas(res))/2; %endpoints
cifG = (sum(midpoints2.*ds)+...
    ds(1)*rnorms(1)+ds(res-1)*rnorms(res))/(4*pi);
%calculate c2 and K
c2 = 1+real(integral);
k = c2+sqrt(1+c2);

%calculate the Cauchy Transform along delOmega
cif = resolvent_norm_integral(A, nr(ind1)) + cifG;

%% plot the numerical range, eigs, and delOmega
figure()
subplot(2,1,1)
plot(nr, '--k'), daspect([1,1,1]), hold on
eigvs = eig(A); 
plot(real(eigvs), imag(eigvs), 'kx');
plot(delOm, '-b', 'LineWidth',2)
%plot e^(At) with k and cif upper bounds
subplot(2,1,2)
% calculate matrix envelope
timestep = 1; iterations = 500;
meA = times_expm(A, timestep, iterations);
tmax = iterations*timestep;
t = 0:timestep:tmax;
semilogy(t,meA, 'DisplayName', 'Matrix Envelope'), hold on
xlim([0,tmax])
plot([0,tmax], [k,k], 'DisplayName', 'K')
plot([0,tmax], [cif, cif], 'DisplayName', 'Cauchy Transform')
legend()