%Natalie Wellen
%10/08/21
%
% Function to estimate the exact value of c1
%input, del_Omega, complex vector, the boundary of the spectral set
%input, sigma_0, complex number, the starting point on del_Omega
%input, sigma_0_prime, complex number, the tangent angle or derivative of
%       del_Omega at sigma_0 in the clock-wise direction
%output, c1, double, the constant c1 associated with the spectral set del_Omega
function c1 = calc_c1(del_Omega, sigma_0, sigma_0_prime)
    %first find where sigma_0 is on del_Omega and set sigma_0 to be the
    %first and last entry 
    start = find(del_Omega == sigma_0, 1);
    del_Om = [del_Omega(start:end), del_Omega(1:start)];
    %calculate the vector of angles between each pair of points on
    %del_Omega
    del_Om = del_Om - sigma_0; %center all the angles at 'zero'
    angles = angle_between(del_Om(2:end-2), del_Om(3:end-1));
    %Add the end_points to the calculation
    a2 = abs(sigma_0_prime - mod(angle(del_Om(end-1)), 2*pi));
    a1 = abs(mod(sigma_0_prime+pi, 2*pi) - mod(angle(del_Om(2)), 2*pi));
    %HOW DOES CHEBFUN CHOOSE THE DIRECTION OF THE DERIVATIVE??
    %From previous examples I think it always points clock-wise, same as
    %how sigma_0_prime is defined in tangent_theta()
    
    %sum up the angles to get c1
    c1 = sum(angles)+a1+a2;
end