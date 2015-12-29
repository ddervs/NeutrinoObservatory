%% Make into 100x100x27
clear all
clc

% The computation time for the location probabilities script scales
% exponentially, so to save time, a 10x10x27 array was made in that script.
% This script takes that data and expands it into a 100x100x27 array using
% a chi^2 fit, as the interaction ength doesn't vary in the x-y directions,
% only in z.

%% Parameters

ndep = 27;
nrow = 10;
ncol = 10;

X = 1:nrow;
Y = 1:ncol;

load('antarctica_10x10_prob.mat') % import data (similar file name for Greenland)

coefficients = cell(27,1);
W = ones(10,10);


%% Chi^2 Fit

for z = 1:ndep
    
    for x = 1:nrow
        for y = 1:ncol
            
            
            Z = normalised_probabilities(:,:,z);
            
            coefficients{z,1} = polyfitweighted2(X,Y,Z,2,W); % Fit x-y plane for each z
            
        end
    end
    
    
end




new_probabilities = zeros(100,100,27);

%% Populate new array

for k = 1:ndep
    
    P = coefficients{k,1};
    
    for i = 1:100
        for j = 1:100
            
            ii = i/10; % rescale
            jj = j/10;
            
            new_probabilities(i,j,k) = P(1)+ii*P(2)+jj*P(3)+(P(4)*ii^2)+(P(5)*ii*jj)+(P(6)*jj^2);
            
        end
    end
end



sperl8 = new_probabilities;

new_probabilities = new_probabilities./sum(sum(sum(new_probabilities)));


%% Plot

figure(1)
myaa
q = slice(sperl8, 1:size(sperl8,1), 1:size(sperl8,2), 1:size(sperl8,3));
alpha('color')
set(q,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
xlabel('width of detector [x100m]')
ylabel('length of detector [x100m]')
zlabel('depth of detector [x100m]')
title('probability of an interaction at any location')