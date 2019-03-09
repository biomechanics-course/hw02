% testfindmin
%
% Unit test for findmin function
% Apply a quadratic function and test whether findmin minimizes it.

% The quadratic 1/2*x'*Q*x + b'*x = 0 has minimum at
% Q*x = b, or x = -Q\b.
Q = [4 -1; -1 2];
b = [1; -1];

% TZ on 03-08-2018: removing the dots from the previous version, such that the outcome is a scalar
quadratictest = @(x) 1/2*x(:)'*Q*x(:) + b'*x(:);

xstar = -Q\b;

% Here is the minimum, using Matlab's fminunc function
xstarfminunc = fminunc(quadratictest, [0;0]);

% Test findmin
xstarfindmin = findmin(quadratictest, [0;0]);

atol = 1e-6; % an absolute tolerance for testing solution
assert(rms(xstarfindmin - xstar) < atol, 'findmin does not meet error tolerance');

% Here is a demonstration of the quadratic function and the minima
[X,Y] = meshgrid(linspace(-5,5),linspace(-5,5));
Z = zeros(size(X));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = quadratictest([X(i,j), Y(i,j)]);
    end
end
clf;
contour(X, Y, Z);
xlabel('x1'); ylabel('x2');
title('Contours of quadratictest'); 
hold on
h = plot(xstar(1), xstar(2), '*', ...
    xstarfminunc(1), xstarfminunc(2), 'o');
legend(h, 'xstar actual', 'xstar findmin');

