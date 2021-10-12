%  User enters matrix A.  This routine computes its field of values
%  and then discretizes it by arclength and takes its derivative.
%  (Note: chebfun stores the boundary of the fov in a clockwise direction;
%  one can reverse the order to make it counter-clockwise, but I don't do that here.)

L = fov(A);   % Compute field of values.
plot(L,'-b'), axis equal, title('Field of values'), shg, pause

%To define the chebfun itself so that it is in counter-clockwise order
S = arcLength(L);   % Measure the length of the bndy.
s = cumsum(abs(diff(L)));
u = inv(s);
u = inv(S-s);    % Reverse the order of L if it is clockwise.
Lnew = newDomain(L,minandmax(u));
W = Lnew(u);     % Now W is bndy of fov discretized by arclength.
Wprime = diff(W);  % Wprime is dW/ds.

%  Now if you want a vector of discrete points on the bndy, do this.
nptsm1 = 49;
ds = S/nptsm1;
Wvec = W([0:nptsm1]*ds);
Wprimevec = Wprime([0:nptsm1]*ds);
%  Now Wvec contains 50 bndy points and Wprimevec contains their derivatives.
hold on; plot(Wvec,'or'), xlabel('Points equally spaced in arc length (red o)'), shg, pause

% Plot the tangent lines.
%for l=1:length(Wprimevec),
%  plot([real(Wvec(l)-Wprimevec(l)),real(Wvec(l)+Wprimevec(l))], ...
%       [imag(Wvec(l)-Wprimevec(l)),imag(Wvec(l)+Wprimevec(l))], '-g');
%end;
%ylabel('Tangent vectors (green)'), shg

%  If you want to reconstruct W from Wvec, just say
Wnew = chebfun(conj(Wvec'),[0,S]);  % I think chebfun wants a column vector and Wvec is a row vector.
plot(Wnew,'-k'), title('Reconstructed fov (black)'), shg, hold off