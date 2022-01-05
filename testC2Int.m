%Script with tests for c2-int




%% 1. Test calc_c2 with A1
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
figure(), plot(del_OmA1), daspect([1,1,1])

%[c2A1] = calc_c2(A1, nrA1, del_OmA1, del_OmA1_prime, 1, 100)
[kA1, c1A1, c2A1] = calc_k(A1, nrA1, del_OmA1, del_OmA1_prime, 1, 100)
%% 2. Test calc_c2 with A1 and two regions removed
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
vertical_index = (real(del_OmA1)> -1/4);
del_OmA1(vertical_index) = -1/4*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = -1i;
vertical_index = (real(del_OmA1)< -1);
del_OmA1(vertical_index) = -1*ones(1,sum(vertical_index)) + 1i*imag(del_OmA1(vertical_index));
del_OmA1_prime(vertical_index) = 1i;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 100)
%% 3. Test calc_c2 with A1 and two regions removed that don't exclude del_Om(1)
A1 = [.5 -1; 1 -1.5];
[nrA1, nrA1_prime] = numerical_range(A1, 100);

% define del_Om for testing
del_OmA1 = nrA1;
del_OmA1_prime = -1i*(del_OmA1 +1/2);
horizontal_index = (imag(del_OmA1)> 0.75);
del_OmA1(horizontal_index) = 0.75*ones(1,sum(horizontal_index))*1i + real(del_OmA1(horizontal_index));
del_OmA1_prime(horizontal_index) = 1;
horizontal_index = (imag(del_OmA1)< -0.75);
del_OmA1(horizontal_index) = -0.75*ones(1,sum(horizontal_index))*1i + real(del_OmA1(horizontal_index));
del_OmA1_prime(horizontal_index) = -1;
figure(), plot(del_OmA1), daspect([1,1,1])

[c2A1] = calc_c2(A1, nrA1, nrA1_prime, del_OmA1, del_OmA1_prime, 1, 100)