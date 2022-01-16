% This function calculates the 2-norm pseudospectral value of epsilon for a
% rectangular region on the complex plane that is then ready to be plotted.
% 
% [epss, wOfPseudo, X, Y] = radius_explore(A, region, resolution)
%  input A, square matrix
%  input, region, [xmin xmax ymin ymax]
%  input, resolution, [xres yres] or res which assumes it is the same along both axes
%  output, epss, real matrix that has dimension xres by yres 
%  output, wOfPseudo, real matrix of the raius of the numerical range for (A-omI)
%  output, X, vector of doubles, real coordinates of grid points
%  output, Y, vector of doubles, imaginary coordinates of grid points
% 
%  Depends on: numerical_range

%Natalie Wellen
%10/26/21

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
            [epss(ii,jj), wOfPseudo(ii,jj)] = r_of_A(A, m, grid(ii,jj)); %num radius
        end
    end
end

%Notice that both outputs are the radius choices for removing a circle
%  centered at grid(ii,jj)