%{
Function to calculate the numerical range of a matrix and output a vector
of finite length

[fov] = numerical_range(A, resolution)
 input, A, the square matrix we are computing the numerical range of
 input, resoltuion, the number of points we are using to estimate the
        numerical range
output fov, vector of complex values depicting the boundary of the
 numerical range
%}

%Natalie Wellen
%10/26/21

function fov = numerical_range(A, resolution)
    %check that A is square
    [n,m] = size(A);
    assert(n==m, 'A must be a square matrix.')
    fov = zeros(1,resolution+1);
    count = 0;
    for j  = [linspace(0, 2*pi, resolution),2*pi]
        count = count+1;
        A_rotated = exp(1i*j)*A;
        B = 1/2*(A_rotated + A_rotated');
        [V,D] = eig(B);
        e = real(diag(D).');
        contributer = max(e);
        %the min real eigenvalue gives the min numerical range rather than max
        index = find(e==contributer);
        fov(count) = V(:,index(1))'*A*V(:,index(1));
    end
end

%Link to algorithm source: https://www.jstor.org/stable/2156587