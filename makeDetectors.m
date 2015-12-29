function [ New_det_loc ] = makeDetectors(num_on_edge, spacing, maxdepth)

%% FUNCTION TO CREATE A HEXAGONALLY TILED ARRAY OF DETECTOR CO-ORDINATES

% This function generates coordinates for a detector array where the
% detectors are arranged in a hexagonal pattern. num_on_edge is the number
% of detectors on each edge, spacing is the distance (in 100m) between them
% and maxdepth (100m) is the depth of the detector.

% Unfortunately this function was not used in the main simulation as a
% square grid pattern showed the most detected events. This script is
% included as a reference however, as it is mentioned in the main text.


New_det_loc = [];

det_on_each_edge = num_on_edge;
distance = spacing;
z = maxdepth;

rows = (det_on_each_edge * 2) - 1;

nascend = 1:rows; % Setting up plot geometry
diffs = abs(nascend - det_on_each_edge);
nrows = rows - diffs;

for i = 1:rows % loop each row
    
    for j = 1:nrows(i) % number of detectors per row
        
        x = ((j-1)+(abs(diffs(i)/2)))*(distance); % place detectors
        y = i*distance;
        New_det_loc = [New_det_loc; [x y z]]; % add to locations
        
    end
    
end


%% OFFSET TO BE CENTRED OVER x = 50, y = 50;

numpoints = length(New_det_loc(:,1));
medindex = median(1:numpoints);
central_point = [New_det_loc(medindex, 1:2) 0];

offset = [];
for k =1:numpoints
    
    offset = [offset; (central_point - [50 50 0])];
    
end

New_det_loc = New_det_loc - offset;

% % PLOT
% 
% X = New_det_loc(:,1);
% Y = New_det_loc(:,2);
% 
% scatter(X, Y);
% 
% axis square
% 


end

