function [ conecoords ] = makeConeTrunc(intstats, cherenkov_angle, detector_depth)

%% makeConeTrunc
% Function takes a row of interaction_stats and Cherenkov angle and 
% generates the coordinates of the corresponding Cherenkov cone. This 
% essentially does the same as makeCone but only takes the top of the cone
% and at higher resolution so as to reduce computation time




xpos = intstats(1); % Take interaction data
ypos = intstats(2);
zpos = intstats(3);
xvel = intstats(5);
yvel = intstats(6);
zvel = intstats(7);

conelength = 15; % Set length of cone

Rrrr = conelength * tan(cherenkov_angle); % Radius of cone

% Cone starts pointing along the z-axis

r=[0:0.05:Rrrr]; % points down length of cone down axis
theta=[0:pi/100:2*pi]; % round the axis
[R,THETA]=meshgrid(r,theta); % make mesh

velmod = sqrt(xvel^2 + yvel^2 + zvel^2);
dir = [(xvel/velmod) (yvel/velmod) (zvel/velmod)];

X=R.*cos(THETA); % make cone points
Y=R.*sin(THETA);
Z=R/tan(cherenkov_angle); 

r = vrrotvec([0 0 1], dir); % rotation vector to shower direction
Rtot = vrrotvec2mat(r); % corresponding rotation matrix

x = [];
y = [];
z = [];

% Put mesh data into vectors to perform rotation transformation 

for i =1:length(X(:,1))
    for j = 1:length(X(1,:))
        
        x = [x X(i,j)];
        y = [y Y(i,j)];
        z = [z Z(i,j)];
        
    end
end

coords = [x;y;z]; % Array of vectors
coords = Rtot*coords; % Rotate

xoffset = xpos * ones(1,length(x)); % offset to interaction location
yoffset = ypos * ones(1,length(y));
zoffset = zpos * ones(1,length(z));

offset = [xoffset; yoffset; zoffset];



conecoords = coords + offset; % Array of x,y,z coordinates of cone

ccdummy = [];

for i = 1:length(conecoords(3,:))

    z = conecoords(3,i); 
    zmax = detector_depth + 0.5;
    zmin = detector_depth - 0.5;
    
    if z <= zmax && z >= zmin  % Only take cone coordinates close to surface
        
        ccdummy = [ccdummy, conecoords(:,i)];
    end
    

end

conecoords = ccdummy;

