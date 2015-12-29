%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  	DETECTOR MONTE-CARLO                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% This script runs the main simulation in a loop repeatedly, which each
% time generates a different detector scenario. This records how many
% neutrino events are observed by the detector, finds an average,
% uncertainty and plots the results.

numofobservedints = [];
numofintstakenplace = [];
simulits = 100000;
rootN = 1/sqrt(simulits);


for asdasf = 1:simulits
    
    
    %% MAIN SIMULATION RUNS UP TO THE POINT BEFORE MULTILATERATION %%
    
    % code omitted to save space
    
    if ~isempty(Detector_Data)
        numofobservedints = [numofobservedints; length(Detector_Data(:,1))];
    else
        numofobservedints = [numofobservedints; 0 ];
        
    end
    
    numofintstakenplace = [numofintstakenplace; Totalinteractions];
    
    clc
    display(['iteration: ', num2str(asdasf)])
    clearvars -except numofobservedints numofintstakenplace rootN simulits
    
    
    
    
end

clc

display(['Expected number of observed interactions/yr: ', num2str(mean(numofobservedints)), ' +- '...
    num2str((std(numofobservedints)*rootN))])

display(['number of interactions that take place/yr: ', num2str(mean(numofintstakenplace)), ' +- ',...
    num2str((std(numofintstakenplace)*rootN))])

hist(numofobservedints, 0:round(max(numofobservedints)))
xlabel('Expected number of observed interactions: ')
ylabel('Frequency')
title(['Greenland: Results of Monte-Carlo simulation, # of iterations: ', num2str(simulits) ])

