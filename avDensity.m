function [ avdensity ] = avDensity(start_depth, end_depth, coeff )

% Function to calculate the average density over a given depth range, using 
% the tanh parametrisation given in the main text.


%% UNITS IN METRES

den = coeff(1);
A = coeff(2);
gamma = coeff(3);
beta = coeff(4);

intAtB = den*start_depth + A*(start_depth+(log(cosh(beta-gamma*start_depth))/gamma));
intAtA = den*end_depth + A*(end_depth+(log(cosh(beta-gamma*end_depth))/gamma));

avdensity = abs((intAtB - intAtA)/(end_depth-start_depth))+1000;

end

