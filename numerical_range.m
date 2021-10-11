%Natalie Wellen

%input A, the matrix we are computing the numerical range of
%input resoltuion, the number of points we are using to estimate the
%  numerical range
%output fov, vector of complex values that "draw" the boundary of the
%  numerical range
function fov = numerical_range(A, resolution)
    fov = [];
    step = 2*pi/(resolution -1);
    for j  = 0:step:2*pi
        A_rotated = exp(1i*j)*A;
        B = 1/2*(A_rotated + A_rotated');
        [V,D] = eig(B);
        e = real(diag(D).');
        contributer = max(e);
        %the min real eigenvalue gives the min numerical range rather than max
        index = find(e==contributer);
        foo = V(:,index(1))'*A*V(:,index(1));
        fov = cat(2, fov, foo);
    end
end

%Link to algorithm source: https://www.jstor.org/stable/2156587