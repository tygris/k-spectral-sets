%Function to define a spectral set Omega by removing disks from W(A), then
% estimate an uppper bound on the value of K for the set Omega.
%
%[k, resNorm, delOm, delOm_prime, c1, c2] = remove_disks(A)
%  input, A, n by n double, the matrix we are defining a spectral set for
%  
%  output, k, double, the K value of the spectral set Omega
%  output, resNorm, double, the value of the integral of the resolvent norm of A
%  output, delOm, complex double, the boundary of Omega in the complex plane 
%  output, delOm_prime, complex double, the derivative of delOm in the
%          counter-clockwise direction
%  output, c1, double, constant used to calculate K
%  output, c2, double, constant used to calculate K

% Depends on:  - numerical_range
%              - nr_gap_fill
%              - cellmat2plot
%              - define_del_Omega
%                  - r_of_A
%                  - delOm_flipper
%              - calc_k_removed_disks
%                  - resolvent_norm_integral
%                  - calc_c1
%                  - calc_c2_curve

%Natalie Wellen
%3/06/23

function [k, resNorm, delOm, delOm_prime, c1, c2] = remove_disks(A)
    %The resolution of W(A). 
    res = 20000;
    %define the starting Omega
    [nr, nr_prime] = numerical_range(A,res);
    [nr, nr_prime] = nr_gap_fill(nr, nr_prime);
    figure()
    plot(nr), daspect([1,1,1])
    %Ask the user for centers of disks to remove from the current state of
    %Omega until they are satisfied
    Y = 'Y'; y = 'Y'; N = 'N'; n = 'N';
    zeta = [];
    r1orr2 = [];
    xs = [];
    more = 'Y';
    delOmCell = {nr};
    delomCell = {zeros(1,length(nr))};
    moveon = 0;
    while more == 'Y'
        zetaNew = input("Where would you like to remove a disk(s)?\n");
        close
        delOm = cellmat2plot(delOmCell,1);
        figure()
        plot(real(delOm), imag(delOm))
        daspect([1,1,1]);
        hold on
        plot(real(zetaNew), imag(zetaNew), 'mo');
        moveon = input('Is this where you would like to remove the disk? (Y/N)\n');
        if moveon == 'Y'
            close
            [delOmNew, delomNew, xsNew, r1orr2New] = define_del_Omega(delOmCell, delomCell, A, zetaNew, res);
            if 'Y' == input('Are you satisfied with the result? (Y/N)\n')
                delOmCell = delOmNew;
                delomCell = delomNew;
                zeta = cat(2, zeta, zetaNew);
                r1orr2 = cat(2, r1orr2, r1orr2New);
                xs = cat(2, xs, xsNew);
            else
                close
                cellmat2plot(delOmCell, delomCell);
            end
        end
        more = input('Would you like to remove more disks? (Y/N)\n');
        if more ~= 'Y'
            if 'N' == input('Are you sure? (Y/N)\n')
                more = input('Would you like to remove more disks? (Y/N)\n');
            end
        end
    end

    close
    %define delOm as a vector
    [delOm, delom] = cellmat2plot(delOmCell, delomCell,1);
    figure()
    plot(nr, '--k'), hold on, axis equal
    plot(delOm,'b')
    %with the results of defining delOm with define_del_Omega(), calculate
    %K and the integral of the resolvent norm
    [k, resNorm, delOm_prime, c1, c2] = calc_k_removed_disks(A, zeta, nr, nr_prime, delOm, delom, xs, r1orr2);
end