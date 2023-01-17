%Function to calculate the spectral set and comparison Cauchy Integral for
% a continuous time dynamical system represented as a matrix
%
%[k,cif,c2] = contDS(A, timestretch)
% input, A, n by n complex double, the continuous time DS
% input (opt), timestretch, double, to shorten or widen the time steps and
%              thus tmax for the second plot
%
% output, k, double, the estimated K value of the spectral set delOm
% output, cif, double, the value of the integral of the resolvent norm
%              along delOmega
% output, c2, double
% output, plot of numerical range overlayed by delOmega and ||exp(At)||
%         bounded by K and the value of cif

%Natalie Wellen
%1/26/22

function [k,cif,c2, cifG, distG, extra_cif] = contDS(A, timestretch)
    if nargin == 1
        timestretch = 1;
    end
    %calculate the numerical range and it's derivative
    [nr] = numerical_range(A,200000);
    %find Gamma1 and it's derivative
    [y1, y2] = nrCutOff(A, 0);
    [c2, cifG] = calc_c2_v(A, imag(y1), imag(y2));
%    mean_cif = cifG/abs(y2-y1);
%    mean_c2 = c2/abs(y2-y1);
    distG = abs(y2-y1);
    %define delOmega and it's derivative
    ind1 = real(nr)<0;
    delOm = cat(2, [y1,y2], nr(ind1));
    %calculate k, c1, and c2
    k = c2 + sqrt(1+c2^2);
    %calculate the Cauchy Transform along delOmega
    extra_cif = cauchyIntFormula(A, nr(ind1));
    cif = extra_cif + cifG;
    
    figure()
    opts= {'LineWidth',2};
    plot(nr, '--k'), daspect([1,1,1]), hold on
    eigvs = eig(A); 
    plot(real(eigvs), imag(eigvs), 'kx');
    plot(delOm, '-b', opts{:})
    title("Spectral Set and Eigenvalues")
    
    %plot the numerical range, eigs, and delOmega
    opts= {'LineWidth',2};
    figure()
    subplot(2,2,2)
    plot(nr, '--k'), daspect([1,1,1]), hold on
    eigvs = eig(A); 
    plot(real(eigvs), imag(eigvs), 'kx');
    plot(delOm, '-b', opts{:})
    title("Spectral Set and Eigenvalues")
    %plot e^(At) with k and cif upper bounds
    
    subplot(2,2,1)
    % calculate matrix envelope
    timestep = 0.1*timestretch; iterations = 500;
    meA = times_expm(A, timestep, iterations);
    tmax = iterations*timestep;
    t = 0:timestep:tmax;
    plot(t,meA, 'DisplayName', 'Matrix Envelope', opts{:}), hold on
    xlim([0,tmax])
    %plot([0,tmax], [k,k], 'DisplayName', 'K', opts{:})
    %plot([0,tmax], [cif, cif], 'DisplayName', 'Cauchy Transform', opts{:})
    %legend()
    
    % plot the resolvent norm and abs(min(eig(A)))
    n = 20001;
    m = length(A);
    Gam1 = linspace(y1, y2, n);
    Gam1_prime = 1i*ones(1,n);
    gammas = zeros(1,n);
    rnorms = zeros(1,n);
    for jj = 1:n
        R = 1/(2*pi)*inv(Gam1(jj)*eye(m) - A); %the matrix s'(sI-A)^-1
        gammas(jj) = -1*(min(eig(R+R'))); 
        rnorms(jj) = norm(R,2);
    end
%    max_rnorm = max(rnorms);
    subplot(2,2,3)
    semilogy(imag(Gam1), rnorms, opts{:})
    ylim([min(rnorms)/5, max(rnorms)*1.5])
    title("Resolvent Norm on Imaginary Axis")
    subplot(2,2,4)
    semilogy(imag(Gam1), gammas, opts{:})
    ylim([min(gammas)/5, max(gammas)*1.5])
    title("|\lambda_{min}| on Imaginary Axis")
end