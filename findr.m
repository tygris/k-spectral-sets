%Function to locate the maximum value of r for a point on a curve with a
%given derivative value. The matrix A must also be given to calculate r
%according to Theorem 2 in "K-Spectral Sets"
%
%[r, piOr2pi] = findr(A, sigma0, sigma0_prime, length, resolution)
% input, A, square matrix we are calculating c2 for
% input, sigma0, complex value, point on del_Omega we are analyzing
% input, sigma0_prime, complex_value, the clockwise direction of the
%        tangent curve to del_Omega at sigma0
% input, length, double, the distance from sigma0 we are searching for a
%        maximum r
% input, resolution, integer, the number of equally spaced points to search
% output, r, double, the length of the radius of the disk tangent to sigma0
%         and with equal derivative (mod pi)
% output, r1orr2, 1 or 2, indicates if 

%Natalie Wellen
%10/25/21

function [r, r1orr2] = findr(A, sigma0, sigma0_prime, length, resolution)
    if ~exist('length', 'var')
        length = 1;
    end
    if ~exist('resolution', 'var')
        resolution = 100;
    end
    search_direction = angle(sigma0_prime) + pi/2; %does it matter if the angle is negative?
    search_points = linspace(sigma0, sigma0 + length*exp(1i*search_direction), resolution);
    
    %set-up for calculating the radii
    [n,m] = size(A);
    if n ~= m
        display("ERROR: A must be a square matrix.")
        return
    end
    r1_vec = zeros(1, resolution);
    r2_vec = zeros(1, resolution);
    for ii = 1:resolution
        A_shift_inv = inv(A - search_points(ii)*eye(n));
        %calculate r1 values (eps pseudo-spectra => min eigv is -1/2*pi*r1)
        r1_vec(ii) = (norm(A_shift_inv))^-1;
        %calculate r2 values (abscissa of the eps p-s => min eigv is -1/pi*r2)
        r2_vec(ii) = 1/(max(abs(numerical_range(A_shift_inv, 100))));
    end
    distance = abs(search_points- sigma0);
    r1_max_index = find(r1_vec >= distance, 1, 'last');
    r2_max_index = find(r2_vec >= distance, 1, 'last');
    r1 = r1_vec(r1_max_index);
    r2 = r2_vec(r2_max_index);
    if 1/(2*pi*r1) >= 1/(pi*r2)
        r1orr2 = 1;
        r = distance(r1_max_index);
    else
        r1orr2 = 2;
        r = distance(r2_max_index);
    end
end