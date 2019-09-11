% cut-off function
function [eta, etap, etapp] = cutoff(delta, l, h)
% Function that computes the cutoff function, and derivatives for the given
% configuration

% Parameters:
%
% INPUT delta, real, given in the .msh file. Upper bound for which the cutoff = 1
% INPUT l, real, given in the .msh file. Lower bound for which the cutoff = 0
% INPUT h, real, meshsize
% OUTPUT eta, function of the transition part
% OUtPUT etap, first derivative of eta
% OUTPUT etapp, second derivative of eta

 %delta = 0.1;
 %l = 0.9;
 %h = 0.001;
 %x = [delta:h:l];
xd = [0: h: delta];
xr = [delta+h: h:l]; 
xl = [l+h: h:1];

xx = unique([xd xr xl]);
%size(x)
L  = length(xx);

a = 6 ./ (l * 0.7 * sqrt(2));
b = ((delta * 6 /l) + 2) ./ (0.7 * sqrt(2));

% Create the continuous function that transitions from 1 to 0 and its
% derivatives
eta   = @(x) 1/2 - (1/2)*erf(a*x-b) ;
etap  = @(x)  -a/ sqrt(pi) * exp ( - x.^2) ;
etapp = @(x) - 2 * a^2 * x / sqrt(pi) .* exp(-x.^2) ;

%% Display
 figure;
 plot(xx, eta(xx), '-xb')
end

