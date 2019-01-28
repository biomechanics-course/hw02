function xstar = findmin(f, x0)
% FINDMIN  Finds the minimum of a function, with vector-valued x0.
%  
% xstar = findmin(f, x0) performs a Newton search to find the 
%   root of the gradient of f, and returns xstar,
%   starting with initial guess x0. The
%   function will typically be expressed with a function handle, e.g. 
%   @f.

% Use the gradient of f as a vector-valued function to find
% the root of. Normally this will be a row vector, but for
% compatibility with findroot, express the gradient as a column
% vector. (The same root is found, whether the function is expressed
% as row or column.)
fgradientv = @(x0) %USE SOME EXISTING FUNCTION(f, x0)';

% The root of this gradient satisfies the
% necessary condition for a minimum.
[xstar, cnvrg] = % ADD SOME CODE HERE

if ~cnvrg % failed to converge
    warning('findmin did not converge on root of gradient');
end

if size(x0,2) > 1   % x0 was given to us as a row vector
    xstar = xstar'; % so return xstar in the same shape
end
    
end % findmin