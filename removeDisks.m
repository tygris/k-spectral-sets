%Function to define a spectral set Omega by removing disks from W(A), then
%estimate an uppper bound on the value of K for the set Omega.
%
%[k, cif, delOm, delOm_prime] = removeDisks(A)
%  input, A, n by n double, the matrix we wish to define a spectral set for
%  
%  output, k, double, the K value of the spectral set Omega
%  output, cif, double, the value of the integral of the resolvent norm of A
%  output, delOm, complex double, the boundary of Omega in the complex plane 
%  output, delOm_prime, complex double, the derivative of delOm in the
%          counter-clockwise direction
%  output, c1, double
%  output, c2, double
%
% Depends on:  - numerical_range
%              - nrGapFill
%              - cellmat2plot
%              - define_del_Omega
%                  - r_of_A
%                  - delOmega_flipper
%              - calc_kRemovedDisk

%Natalie Wellen
%1/18/23

function [k, cif, delOm, delOm_prime, c1, c2] = removeDisks(A)
    %The resolution of the W(A). Half this value is used for the
    %disks removed from W(A).
    res_num = 2000;
    %define the starting Omega
    [nr, nr_prime] = numerical_range(A,res_num);
    [nr, nr_prime] = nrGapFill(nr, nr_prime);
    figure()
    plot(nr), daspect([1,1,1])
    %Ask the user for centers of disks to remove from the current state of
    %Omega until they are satisfied
    Y = 'Y'; y = 'Y'; N = 'N'; n = 'N';
    om = [];
    r1orr2 = [];
    xs = [];
    more = 'Y';
    del_Om = {nr};
    del_om = {zeros(1,length(nr))};
    moveon = 0;
    while more == 'Y'
        om_new = input("Where would you like to remove a disk(s)?\n");
        close
        del_Om_vec = cellmat2plot(del_Om,1);
        figure()
        plot(real(del_Om_vec), imag(del_Om_vec))
        daspect([1,1,1]);
        hold on
        plot(real(om_new), imag(om_new), 'mo');
        moveon = input('Is this where you would like to remove the disk? (Y/N)\n');
        if moveon == 'Y'
            close
            [del_Om_new, del_om_new, xs_new, r1orr2_new] = define_del_Omega(del_Om, del_om, A, om_new, res_num);
            if 'Y' == input('Are you satisfied with the result? (Y/N)\n')
                del_Om = del_Om_new;
                del_om = del_om_new;
                om = cat(2, om, om_new);
                r1orr2 = cat(2, r1orr2, r1orr2_new);
                xs = cat(2, xs, xs_new);
            else
                close
                cellmat2plot(del_Om, del_om);
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
    [delOm, delomvec] = cellmat2plot(del_Om, del_om,1);
    figure()
    plot(nr, '--k'), hold on, axis equal
    plot(delOm,'b')
    %with the results of defining delOm with define_del_Omega(), calculate
    %K and the integral of the resolvent norm
    [k, cif, delOm_prime, c1, c2] = calc_kRemovedDisk(A, om, nr, nr_prime, delOm, delomvec, xs, r1orr2);
end