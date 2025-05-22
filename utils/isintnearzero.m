function f = isintnearzero(m)
% function f = isintnearzero(m)
%
% <m> is a matrix
%
% return a logical matrix the same size as <m>.
% an element is 1 iff it is a float, finite, equal to an integer with some 
% tolerance to some imprecision (default tolerance = 1.0^-10):
% 
%   f = isfloat(m) & isfinite(m) & nearZero(m-round(m));
%
% example:
% isequal(isintnearzero([1 1.5 NaN Inf]),[1 0 0 0])
%
% history:
% 03/2025: take original isint.m from knkutils and add nearZero check.
% 05/2025: change name from isint.m to isintnearzero.m

% do it
f = isfloat(m) & isfinite(m) & nearZero(m-round(m));
