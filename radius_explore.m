%This function calculates the 2-norm pseudospectral value of epsilon for a
% rectangular region on the complex plane
%
%[epss, wOfPseudo, X, Y] = radius_explore(A, region, resolution)
% input A, square matrix
% input, region, [xmin xmax ymin ymax]
% input, resolution, [xres yres] or res which assumes it is the same along both axes
% output, epss, real matrix that has dimension xres by yres 
% output, wOfPseudo, real matrix of the raius of the numerical range for (A-omI)
% output, X, the grid spacing on the real axis (for plotting)
% output, Y, the grid spacing on the imaginary axis (for plotting)

%Natalie Wellen
%10/25/21

function [epss, wOfPseudo, X, Y] = radius_explore(A, region, resolution)
    %check the dimension of A
    [m,n] = size(A);
    if m~= n
        disp("A needs to be a square matrix")
        return
    end
    
    %create the grid of values to calc eps for
    numRes = size(resolution,2);
    if numRes == 1
        resolution = [resolution resolution];
    end
    xvals = linspace(region(1), region(2), resolution(1));
    yvals = linspace(region(3), region(4), resolution(2));
    [X,Y] = meshgrid(xvals, yvals);
    grid = X+1i*Y;
    
    %define the output variable to save values in
    epss = zeros(resolution(1), resolution(2));
    wOfPseudo = zeros(resolution(1), resolution(2));
    
    %move through the meshgrid and calculate epsilon and numerical readius
    for ii = 1:resolution(1)
        for jj = 1:resolution(2)
            A_shift_inv = inv(A-grid(ii,jj)*eye(m));
            epss(ii,jj) = norm(A_shift_inv); %epsilon
            wOfPseudo(ii,jj) = max(abs(numerical_range(A_shift_inv,100))); %num radius
        end
    end
    epss = epss.^-1; %radius with lam_min bounded by -R/2pi
    wOfPseudo = wOfPseudo.^-1; %radius with lam_min bounded by -R/pi
end

%Notice that both outputs are the radius choices for removing a circle
%  centered at grid(ii,jj)