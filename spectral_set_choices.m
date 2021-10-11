%Natalie Wellen
%10/08/21
%
%This function creates two contour plots of the potnential radii of circles
%that can be removed from the numerical range with bound -R/2pi or -R/pi on
%the minimum eigenvalues of the operator mu
%
%Input A, square matrix that we are interested in examining the spectral
%  set of 
%Input numRange, complex vector, boundary of the numerical range of A
%Input resolution, integer of 1 by 2 vector of integers, the number of 
% points in the x and y direction used to define the grid for calculating
% potential radii 
%Output epss, matrix of points of the psuedospectral value of A at x+iy
%(-R/2pi)
%Output wOfPseudo, matrix of the numerical radius of inv(A - x+iyI)
%(-R/pi)
%Output, contour plot of epss with the numerical range outlined
%Output, contour plot of wOfPseudo with the numerical range outlined
function [epss, wOfPseudo] = spectral_set_choices(A, numRange, resolution)
    %set boundaries of grid to calculate radii for
    xmin = min(real(numRange));
    xmax = max(real(numRange));
    ymin = min(imag(numRange));
    ymax = max(imag(numRange));
    
    %calculate the potential radii
    [epss,wOfPseudo, X, Y] = radius_explore(A, [xmin, xmax, ymin, ymax], resolution);
    
    %create the plot for epss
    figure()
    [M1,c1] = contour(X,Y,epss,15,'ShowText', 'on','DisplayName', 'Radii');
    daspect([1,1,1])
    c1.LineWidth = 2;
    colorbar();
    title('Radius of circle to be removed with \lambda_{min} \geq -R/(2\pi)')
    %make y-ticks imaginary?
    hold on 
    plot(numRange, '-m', 'LineWidth', 1.5, 'DisplayName', 'Numerical Range')
    legend('location', 'southoutside','NumColumns',2)
    hold off
    
    %create the plot for wOfPseudo
    figure()
    [M2,c2] = contour(X,Y,wOfPseudo,15,'ShowText', 'on','DisplayName', 'Radii');
    c2.LineWidth = 2;
    colorbar();
    title('Radius of circle to be removed with \lambda_{min} \geq -R/\pi')
    daspect([1,1,1])
    %make y-ticks imaginary?
    hold on 
    plot(numRange, '-m', 'LineWidth', 1.5, 'DisplayName', 'Numerical Range')
    legend('location', 'southoutside','NumColumns',2)
    hold off
end