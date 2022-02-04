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

%uses inpolygon

%Natalie Wellen
%2/02/22

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
    [Gam1, ~] = nrDiskOff(nr,1); %I don't use the points nr intersects D 
    Gam1_prime = 1i*Gam1; %counter-clockwise derivative of the unit disk
    %calculate c2 and k
    [c2, cifG] = calc_c2_curve(A, Gam1, Gam1_prime);
    k = c2 + sqrt(1+c2^2);
    %calculate the Cauchy Transform along delOmega
    cif = cauchyIntFormula(A, nr(ind1)) + cifG;
    
%     %plot the numerical range, eigs, and delOmega
%     opts= {'LineWidth',2};
%     figure()
%     subplot(2,2,2)
%     plot(nr, '--k'), daspect([1,1,1]), hold on
%     evs = eig(A); 
%     plot(real(evs), imag(evs), 'kx');
%     plot(Gam1, '-b', opts{:}), plot(nr(ind1), '-b', opts{:})
%     title("Spectral Set and Eigenvalues")
%     
%     %plot ||A^k|| with K and cif upper bounds
%     subplot(2,2,1)
%     % calculate matrix envelope
%     iterations = 500*time_stretch;
%     meA = times_expm(A, timestep, iterations);
%     tmax = iterations*timestep;
%     t = 0:timestep:tmax;
%     plot(t,meA, 'DisplayName', 'Matrix Envelope', opts{:}), hold on
%     xlim([0,tmax])
%     plot([0,tmax], [k,k], 'DisplayName', 'K', opts{:})
%     plot([0,tmax], [cif, cif], 'DisplayName', 'Cauchy Transform', opts{:})
%     %legend()
%     
%     % plot the resolvent norm and abs(min(eig(A)))
%     n = 20001;
%     m = length(A);
%     Gam1 = linspace(y1, y2, n);
%     Gam1_prime = 1i*ones(1,n);
%     gammas = zeros(1,n);
%     rnorms = zeros(1,n);
%     for jj = 1:n
%         R = 1/(2*pi)*inv(Gam1(jj)*eye(m) - A); %the matrix s'(sI-A)^-1
%         gammas(jj) = -1*(min(eig(R+R'))); 
%         rnorms(jj) = norm(R,2);
%     end
%     subplot(2,2,3)
%     semilogy(imag(Gam1), rnorms, opts{:})
%     ylim([min(rnorms)/5, max(rnorms)*1.5])
%     title("Resolvent Norm on Imaginary Axis")
%     subplot(2,2,4)
%     semilogy(imag(Gam1), gammas, opts{:})
%     ylim([min(gammas)/5, max(gammas)*1.5])
%     title("|lambdamin| on Imaginary Axis")
end








