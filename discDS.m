%Function to calculate the spectral set and comparison Cauchy Integral for
% a discrete time dynamical system represented as a matrix
%
%[k,cif,c2] = discDS(A, timestretch)
% input, A, n by n complex double, the discrete time DS
% input (opt), timestretch, double, to shorten or widen the time steps and
%              thus tmax for the second plot
%
% output, k, double, the estimated K value of the spectral set delOm
% output, cif, double, the value of the integral of the resolvent norm
%              along delOmega
% output, c2, double
% output, plot of numerical range overlayed by delOmega and ||exp(At)||
%         bounded by K and the value of cif
%
%Depends on: -numerical_range
%            -nrDiskOff
%            -calc_c2_d
%            -cauchyIntFormula

%Natalie Wellen
%2/03/22

function [k,cif,c2] = discDS(A, timestretch)
    if nargin == 1
        timestretch = 1;
    end
    %calculate the numerical range 
    nres = 200000;
    [nr] = numerical_range(A,nres);
    %Omega is the intersection of the unit disk and W(A)
    ind1 = abs(nr)<=1; %indices of the numerical range still part of delOm
    %find Gamma1 
    [Gam1, as] = nrDiskOff(nr,1); %I don't use the points nr intersects D 
    if isempty(as)
        as = [0, 2*pi];
    end
    %Gam1_prime = 1i*Gam1; %counter-clockwise derivative of the unit disk
    %calculate c2 and k
    m = length(as);
    c2 = 0; cifG = 0;
    for jj = 1:2:m-1
        [c2hold, cifGhold] = calc_c2_d(A, as(jj), as(jj+1));
        c2 = c2+c2hold; cifG = cifG+cifGhold;
    end
    k = c2 + sqrt(1+c2^2);
    %calculate the Cauchy Transform along delOmega
    if as == [0, 2*pi]
        cif = cifG;
    else
        cif = cauchyIntFormula(A, nr(ind1)) + cifG;
    end
    
    %plot the numerical range, eigs, and 
    opts= {'LineWidth',2};
    figure()
    subplot(2,2,2)
    plot(nr, '--k'), daspect([1,1,1]), hold on
    evs = eig(A); 
    plot(real(evs), imag(evs), 'kx');
    %plot delOmega
    nr_plot = nr;
    nr_plot(~ind1) = NaN;
    plot(Gam1, '-b', opts{:}), plot(nr_plot, '-b', opts{:})
    title("Spectral Set and Eigenvalues")
    
    %plot ||A^k|| with K and cif upper bounds
    subplot(2,2,1)
    % calculate matrix envelope
    iterations = 100*timestretch;
    Aks = normAk(A, iterations);
    t = 1:iterations;
    plot(t,Aks, 'DisplayName', 'norm(A^k)', opts{:}), hold on
    xlim([0,iterations])
    plot([0,iterations], [k,k], 'DisplayName', 'K', opts{:})
    plot([0,iterations], [cif, cif], 'DisplayName', 'Cauchy Transform', opts{:})
    %legend()
    
    % plot the resolvent norm and abs(min(eig(A)))
    angley = @(y) (y+2*pi).*(y<0)+y.*(y>=0); %convert to [0, 2pi] domain
    n = 20001;
    m = length(A);
    Gam1 = exp(1i*(linspace(as(1), as(2), n)));
    Gam1_prime = 1i.*Gam1;
    gammas = zeros(1,n);
    rnorms = zeros(1,n);
    for jj = 1:n
        R = 1/(2*pi)*inv(Gam1(jj)*eye(m) - A); %the matrix s'(sI-A)^-1
        gammas(jj) = -1*(min(eig(R+R'))); 
        rnorms(jj) = norm(R,2);
    end
    subplot(2,2,3)
    semilogy(angley(angle(Gam1)), rnorms, opts{:})
    ylim([min(rnorms)/5, max(rnorms)*1.5])
    title("Resolvent Norm on the Unit Disk")
    subplot(2,2,4)
    semilogy(angley(angle(Gam1)), gammas, opts{:})
    ylim([min(gammas)/5, max(gammas)*1.5])
    title("|\lambda_{min}| on the Unit Disk")
end

function [Aks] = normAk(A, iterations)
    %define output vector
    Aks = zeros(1,iterations);
    Ak = A;
    Aks(1) = norm(Ak, 2);
    for jj = 2:iterations
        Ak = Ak*A;
        Aks(jj) = norm(Ak, 2);
    end
end
