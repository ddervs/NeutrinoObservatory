
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%
%%                             MAIN SIMULATION                           %%
%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All code in this script and others following written by Danial Dervovic 
% and Mark Longhurst, apart from polyfitweighted2.m, which is taken from
% http://www.mathworks.co.uk/matlabcentral/fileexchange/13719-2d-weighted
% -polynomial-fitting-and-evaluation/content/polyfitweighted2/polyfitweig
% hted2.m

% All simulation parameters for Greenland with the Antarctica parameters
% given in comments next to the parameters which change.

clear all;
clc;
close all;

% This is the main simulation in which the activity of a radio neutrino
% observatory is simulated. It consists of 3 parts, Event Generation,
% Detection and Event Reconstruction. The first part simulates the physics
% of how neutrinos interact within the detection medium (ice). The second
% part models how the resulting Askaryan radiation from the interactions is
% measured by the detectors. The third part then uses this data to
% reconstruct the interaction vertex, using only the data that has been
% 'measured'.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            EVENT GENERATION                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the code uses Monte-Carlo methods to calculate the number of
% expected interactions in the detector in one year. Then using this
% number, it generates positions and velocities of these interactions,
% again using Monte-Carlo methods. Monte-Carlo methods are also used to
% generate Energies for the incoming neutrinos.

% The probabilities used to find the positions and velocities come are
% generated in the 'probMatrixwDensity.m' script. These values are then
% used in the next section to generate Cherenkov cones which can be deteted
% by the detector array.

%% Constants

n1 = 100000; %number of iterations
% avogadro's number
Na = 6.0221413E23;
% molecular mass of ice
A = 18.01528;
% % preallocate space

% % Interaction Length of a neutrino

load('100x100x27greenland.mat') % Matrix of probabilities of interaction locations
xsec = 5.3E-32;%cm2 interaction cross section
u = 1.66E-24; %g

detector_width = 10;
detector_length = 10;
detector_depth = 27; % times 100m % FOR ANTARCTICA = 28

%% Density profile

den = 25.75; % Coefficients of density profile % FOR ANTARCTICA den = 300
gam = 0.022;                                    % gam = 0.04
AA = 1.1;                                       % AA = 325
bet = 75;                                       % bet = 1

denscoeff = [den AA gam bet];

rho = avDensity(0,100*detector_depth,denscoeff); %gcm-3 % ANTARCTICA use avRefIndex as no added coefficient
l = u/(xsec*(rho/1000)); %interaction length in cm


%% Probability of interaction location
sperl = new_probabilities;
S = sperl;

figure(1) % Plot location likelihoods
q = slice(sperl, 1:size(sperl,1), 1:size(sperl,2), 1:size(sperl,3));
alpha('color')
set(q,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
xlabel('width of detector [x100m]')
ylabel('length of detector [x100m]')
zlabel('depth of detector [x100m]')
title('probability of an interaction at any location')


hhbar = waitbar(0,'progress');
for jbb = 1:1; % Loop needed to show waitbar
    
    %% Variables
    
    % ice density g/cm3 varying with temp
    Rho = rho/1000 + (5.4E-3.*rand([1,n1]));
    
    % ice depth (cm) randomly +/-50cm
    d = detector_depth*10000+((50-(100.*rand([1,n1]))));
    
    % cosmic neutrino energy (eV)
    elow = 1E15;
    eup = 1E20;
    epower = round(15+(19-15).*rand(1,n1));
    est = (1+(10-1).*rand(1,n1));
    eV =sort(est.*(10.^epower)); % Vector of randomly generated energies
    
    
    %% Fluxes - initally all (100km-2 y-1)
    
    % electron neutrino/antineutrino from atmosphere
    fluxe = (eV./1E9).*(24*3600*365*1E10).*(12.57*7.8E-11).*(eV/1E12).^-3.6;
    % muon neutrino/antineutrino from atmosphere
    fluxm = (eV./1E9).*(24*3600*365*1E10).*(12.57*1.05E-10).*((eV./1E12).^-3.7);
    % AGN-NMB neutrinos
    fluxagn = (eV./1E9).*(24*3600*365*1E10).*(12.57*10^-6).*((eV./1E12).^-3.5);
    %Neutrinos from GZK effect
    alp = -2.0;
    Emax = 10^20.5;%(GeV)
    fluxgzk = (24*3600*365*1E10).*0.5.*((eV.^alp).*exp(-eV./Emax));
    %total flux
    Tflux = fluxe+fluxm+fluxagn+fluxgzk;
    %gives around 9 interactions in the detector per year. to be expected if
    %icecube has detected 1-3 a year and is x100 smaller and only looking at
    %energyies of (1E15-1E16)
    
    %% Interaction Cross-sections
    
    i = 1:n1;
    % interaction cross-section (sqcm) for Charged current reactions
    xseccc(i) = 2.69E-36.*((eV./1E9).^0.402); %no neutral current interactions above 1E12eV
    
    % assume averaged interaction cross-section (sqcm)
    xsectot(i) = (xseccc);
    
    %% Number of interactions
    
    % unit check(cm2*(100km-2)-2*g*g-1*cm-3*cm*y-1)
    Ninte = xsectot.*fluxe.*((Na/A).*Rho.*d);
    Nintm = xsectot.*fluxm.*((Na/A).*Rho.*d);
    Nintagn = xsectot.*fluxagn.*((Na/A).*Rho.*d);
    Nintgzk = xsectot.*fluxgzk.*((Na/A).*Rho.*d);
    Ninttot = xsectot.*Tflux.*((Na/A).*Rho.*d);
    
    Tinte = mean(Ninte);
    Tintm = mean(Nintm);
    Tintagn = mean(Nintagn);
    Tintgzk = mean(Nintgzk);
    
    display('Interactions in the detector per year')
    Totalinteractions = mean(Ninttot)
    
    numintround = round(Totalinteractions);
    
    
    n2 = 100; %input('detector length/width - ');
    m = detector_depth; %input('detector depth - ');
    Det = zeros(n2,n2,m); %decimeters
    
    
    
    for io = 1:1
        for tet = 1:(round(Totalinteractions))%# possible numbers
            aa = 1:n2;
            bb = 1:m;
            w1 = S(:,1,1)'; % Weightings for x and y probabilities
            w2 = permute(S(1,1,:),[3 2 1])'; % Weightings for z probabilities
            N = 1;
            
            i = round(aa( sum( bsxfun(@ge, rand(N,1), cumsum(w1./sum(w1))), 2) + 1 ));
            if i==1
                i = 2;
            else
                i=i;
            end
            if i==100
                i = 99;
            else
                i=i;
            end
            j = round(aa( sum( bsxfun(@ge, rand(N,1), cumsum(w1./sum(w1))), 2) + 1 ));
            if j==1
                j = 2;
            else
                j=j;
            end
            if j==100
                j = 99;
            else j=j;
            end
            k = round(bb( sum( bsxfun(@ge, rand(N,1), cumsum(w2./sum(w2))), 2) + 1 ));
            if k==1
                k = 2;
            else
                k=k;
            end
            if k==27
                k = 26;
            else
                k=k;
            end
            Det(sub2ind(size(Det),i,j,k)) = 1;
            pi(tet) = (i);
            pj(tet) = (j);
            pk(tet) = (k);
        end
    end
    position_of_interactions = [pi',pj',pk']; % positions of interactions
    
    %% Plot of Neutrino Energy vs Flux
    
    figure(3)
    
    bf = polyfit(log(eV), log(Tflux), 1);
    % compute fit in linear domain
    Tflux_hat = exp(bf(1) * log(eV) + bf(2));
    % make log log plot
    loglog(eV, Tflux_hat,'.', eV, Tflux);
    label = ['log(Tflux) = ' num2str(bf(1)) 'log(eV) + ' num2str(bf(2))];
    legend('data', label);
    % hold off;
    
    [best_f] = 10.^(bf(1,1)*log(sort(Tflux))+bf(1,2));%will give weightings to energy of interacting neutrino
    
    
    
    for tet = 1:(round(Totalinteractions))% # possible numbers
        ce = sort(eV);
        w3 = sort(best_f,'descend'); %
        N = 1;              %# how many numbers to generate
        
        ei(tet) = ce( sum( bsxfun(@ge, rand(N,1), cumsum(w3./sum(w3))), 2) + 1 ); % This vector is the energies of the interactions
    end
    
    
    
    
    %% Velocity from relativistic energy
    
    ce = 3E8;
    m0c2 = 1.6E-19; %ev/c^2
    E = ei;
    % E  = (m0*c^2)/(sqrt(1-beta))
    % where beta = (v/c)^2
    
    beta = 1-(((m0c2)./E).^2);
    v = (sqrt(beta)*ce)';
    
    interaction_stat = ([position_of_interactions,v]);
    
    %% Likely direction
    
    for qq = 1:length(position_of_interactions)
        ii = position_of_interactions(qq,1);
        iup = n2;
        a = [1:(ii-1),(ii+1):iup]; % Possible values
        wa = [(ii-1):1,1:(iup-ii)]; % Weighting
        ip(qq) = (a( sum( bsxfun(@ge, rand(1,1), cumsum(wa./sum(wa))), 2) + 1 )); % Where flux is from
        
    end
    
    for qq = 1:length(position_of_interactions)
        jj = position_of_interactions(qq,2);
        jup = n2;
        b = [1:(jj-1),(jj+1):jup]; % Possible values
        wb = [(jj-1):1,1:(jup-jj)]; % Weighting
        jp(qq) = (b( sum( bsxfun(@ge, rand(1,1), cumsum(wb./sum(wb))), 2) + 1 )); % Where flux is from
        
        
    end
    
    for qq = 1:length(position_of_interactions)
        kk = position_of_interactions(qq,3);
        iup = m;
        c = [1:(kk-1),(kk+1):iup]; % Possible values
        wc = [(kk-1):1,1:(iup-kk)]; % Weighting
        kp(qq) = (c( sum( bsxfun(@ge, rand(1,1), cumsum(wc./sum(wc))), 2) + 1 )); % Where flux is from
        
        
    end
    
    %% Monte-Carlo generated interaction data
    
    % The interaction_stats array holds the data for each interaction
    % generated, i.e. its x,y,z (100m) position within the detector media, its
    % speed (m/s) and its velocity direction (un-normalised, in 100m/s^-1).
    % The next part of the code deals with how the detectors react to the
    % E-field generated by the Askaryan effect of the resulting electron
    % shower.
    
    stat = ([ip',jp',kp']);
    interaction_stats = ([interaction_stat,stat]);
    
    
    
    %% Plot interaction locations and paths through ice to get there
    
    figure(4)
    xlabel('width of detector [x100m]')
    ylabel('length of detector [x100m]')
    zlabel('depth of detector [x100m]')
    title('Location and path of Interactions + radiation')
    for ppo = 1:round(Totalinteractions)
        grid on
        hold on
        stats_plotable = [interaction_stats(ppo,1:3);interaction_stats(ppo,5:7)]';
        plot3([stats_plotable(2,1),stats_plotable(2,2)],[stats_plotable(1,1),stats_plotable(1,2)],[stats_plotable(3,1),stats_plotable(3,2)],'--')
        
        plot3(stats_plotable(2,1),stats_plotable(1,1),stats_plotable(3,1),'o')
        
    end
    hold off
    axis([0 100 0 100 0 27])
    daspect([1 1 1])
    waitbar(jbb/10)
end
delete(hhbar)

interaction_stats(:,5:7) = interaction_stats(:,5:7)-position_of_interactions;

% for i = 1:length(interaction_stats(:,1))
%
%    interaction_stats(i,5:7) = interaction_stats(i,5:7) / norm(interaction_stats(i,5:7) );
%
% end

interaction_stats
display('        ------------------------------------------------------------------------------')
display('           x           y           z     v(m/s)            i           j           k ')

intstats = []; % Matrix to be filled with interaction information only for detected interactions





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               DETECTION                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This section of the code proceeds by taking the Monte-Carlo generated
% interactions, and finding which of the interactions would be detected.
% The code works via a process of elimination, eventually returning a cell
% array, Detector_Data. Each row corresponds to one of the DETECTED
% interactions. Each element of the first column of the array
% corresponds to one detector which has received a signal from the neutrino
% interaction. Along each row, (in order) is: the x,y,z co-ordinates of the
% detector, the distance from the interaction, the time the signal takes to
% get to the detector, the time difference from the detector's lowest time
% recorded, the angular deviation from the Cherenkov cone, the Magnitude
% of the E-Field observed at the detector, the x,y,z components of
% polarisattion of the E-field. The second element of each column of the
% cell array corresponds to an event ID (useful for checking the correct
% event was being measured).

% Note that in the final part of the simulation, event reconstruction, only
% the Time difference (TDOA) data and E-Fields are used to reconstruct the
% position, velocity and energy of the interaction, along with the detector
% positions; as this is the data that would be gathered in a real detector.

load('NEW_det_loc.mat') % Array of detector x,y,z co-ordinates

% Detection sphere with r = detection_radius around each detector. If the
% interaction lies within the detection radius then the detector MIGHT
% fire. This depends on how far away from the Cherenkov the detector lies.


%% Detection radius

% Use the characteristics of ARA's radio attenae to give a detection radius

operating_frequency = 500; %MHz
bandwidth = 600;           %MHz
noise_pow_spec_density = 3.981e-21; % (converted to W) given as -174dBm in http://arxiv.org/pdf/1105.2854v2.pdf

white_noise_temp = noise_pow_spec_density/(1.38*10^(-23));


sn_ratio = 5; % demanded signal/noise ratio
f = 2; % from paper

fac = 5 * f * sqrt(sn_ratio*white_noise_temp)*(1+(1.6*10^(-6)*operating_frequency^2)/sqrt(bandwidth));

detrad = 10 * 1e15 * 1e-12 / fac;
detection_radius = detrad;              % A radius value



%%  Make detection spheres
figure(5)
hold on
daspect ([1 1 1])
axis([0 100 0 100 0 detector_depth])
[x,y,z] = sphere;      %# Makes a 21-by-21 point sphere
x = x(1:11,:);       %# Keep top 11 x points
y = y(1:11,:);       %# Keep top 11 y points
z = z(1:11,:);       %# Keep top 11 z points

%% Plot Detection Spheres

Detector_location = New_det_loc;
for j = 1:size(Detector_location)
    dsp = mesh((detection_radius.*x)+(Detector_location(j,1)),(detection_radius.*y)+(Detector_location(j,2)),(detection_radius.*z)+(Detector_location(j,3)));  %# Plot the surface
    
end

for k = 1:numintround
    scatter3(position_of_interactions(k,1),position_of_interactions(k,2),position_of_interactions(k,3),'o','fill','k','LineWidth',6)
end


%% Detector Data - put into a cell for easier looping

Detector_Data = cell(numintround,2);


% This part selects for detectors which are within the detection radius

for k = 1:numintround % Loop through each interaction
    
    Det_location = [];
    for j = 1:size(Detector_location(:,1)) % loop through detectors
        
        distfromdet = sqrt(((position_of_interactions(k,1)-Detector_location(j,1)).^2)+...
            ((position_of_interactions(k,2)-Detector_location(j,2)).^2)+...
            ((position_of_interactions(k,3)-Detector_location(j,3)).^2));
        
        
        if detection_radius >= distfromdet % then interaction is detected
            
            detloc = [Detector_location(j,1) Detector_location(j,2) Detector_location(j,3) distfromdet];
            Det_location = [Det_location; detloc];
            
            
        end
        
        
    end
    
    % Remove interaction if detected at less than 6 detectors (DTOA only
    % works for 5+ detectors)
    % This is done at a later stage as well but shortens computation
    % time to do here also
    
    if sum(size(Det_location)) == 0 % check if no detections
        
    elseif length(Det_location(:,1)) < 6 % check if fewer than 6
        
    else
        Detector_Data{k,1} = Det_location;
        Detector_Data{k,2} = k;
        intstats = [intstats; [interaction_stats(k,:) k]];
    end
    
    
    
    
    
    
end


% Remove empty interactions

emptyCells = cellfun('isempty', Detector_Data);
if sum(emptyCells) > 0
    Detector_Data(emptyCells) = [];
    len = length(Detector_Data)/2;
    Detector_Data = reshape(Detector_Data, round(len), 2);
    newnumint = length(Detector_Data(:,1));
else
    display('no detectable interaction');
    return
end

for i = 1:newnumint
    
    DData = Detector_Data{i,1};
    ddatasize = length(DData(:,1));
    DData = [DData, zeros(ddatasize,2)]; % Empty slots for DTOA information, added here
    % due to way code was originally
    % structured
    
    
    Detector_Data{i,1} = DData;
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% E-Fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the Detection section deals with the transmission of
% E-fields and the Cherenkov effect, which depend on the attenuation length
% and refractive index. These quantities vary throughout the ice. To
% account for this, these quantities are averaged over the paths of
% radiation through the ice to the detectors of interest.


%% Refractive Index

% Parameters for variation of refractive index in the ice

% Make arrays of average refractive index, Cherenkov angle and speed of
% transmission for each interaction

% Parameters

refind0 = 1.3; % FOR ANTARCTICA: refind0 = 1.375
A = 0.225;                     % A = 0.187
gamma = 0.037;                  % gamma = 0.05
beta = 1.1;                     % beta = 1.1

ncoeff = [refind0 A gamma beta];
average_ref_index = [];
avcangle = [];
cspeed = [];


for i = 1:newnumint
    
    dep = intstats(i,3)*100;
    n = avRefIndex(dep, 0, ncoeff);
    average_ref_index(i) = n;
    avcangle(i) = acos(1/n);
    cspeed(i) = 300000000/n;
    
end

% refractive_index_of_ice = 1.4;
% cherenkov_angle = acos(1/refractive_index_of_ice);

%% Attenuation length

% Attenuation length given as a polynomial, integrate to find average
% attenuation length for each path from the interaction to the detectors

x6 = -4E-15; % FOR ANTARCTICA x6 = 0.2012
x5 = 3E-12;  %                x5 = -2.063
x4 = 6E-10;  %                x4 = 8.1838
x3 = 1E-06;  %                x3 = 15.627
x2 = 0.0004; %                x2 = 15.553
x1 = 0.0233; %                x1 = 6.5676
x0 = - 1.4911;%               x0 = 2.5307

attencoeff = [x0 x1 x2 x3 x4 x5 x6];


%% Plot Cherenkov Cones

figure(6)
hold on
daspect ([1 1 1])
axis([0 100 0 100 0 detector_depth])

for i = 1:newnumint
    
    
    cone = makeCone(intstats(i,:),avcangle(i));
    x = cone(1,:);
    y = cone(2,:);
    z = cone(3,:);
    scatter3(x,y,z,0.5,'.');
    
end

%% Find E-Field at Cherenkov Angle

nu_0 = 500;
nu = operating_frequency;
nu_scaled = nu/nu_0;
showerE = ei;
e_cherenkov_times_R = 1.1*10^(-7).*(showerE./(1*10^12))*nu_scaled*(1/(1+0.4*nu_scaled^2));

%% Find Ring Centre

ring_centre = zeros(newnumint,3);
for i = 1:newnumint
    
    interaction_point = intstats(i,1:3);
    interaction_velocity = intstats(i,5:7);
    lambda = (detector_depth - interaction_point(3))/interaction_velocity(3);
    
    if lambda >= 0 % Upward moving neutrinos
        rx = interaction_point(1) + lambda*interaction_velocity(1);
        ry = interaction_point(2) + lambda*interaction_velocity(2);
        ring_centre(i,:) = [rx ry detector_depth];
    else % Downward moving neutrinos
        lambda = (1 - interaction_point(3))/interaction_velocity(3);
        rx = interaction_point(1) + lambda*interaction_velocity(1);
        ry = interaction_point(2) + lambda*interaction_velocity(2);
        ring_centre(i,:) = [rx ry 1];
        
        % Check if ring doesnt hit surface
        vec2check = ring_centre(i,:) - interaction_point;
        % Dot product with z direction
        dp = dot([0 0 1], vec2check);
        angle = acos(dp)/norm(vec2check)-(3.141592/2);
        if angle > avcangle(i)
            ring_centre(i,:) = NaN;
        end
    end
    
    
    
    
    
end

%% Remove interactions where cone doesnt hit surface

intstatsdummy = [];
average_ref_index_dummy = [];
avcangle_dummy = [];
cspeed_dummy = [];
for i = 1:newnumint
    
    
    
    if isnan(ring_centre(i,1))
        
        Detector_Data(i,:) = [];
        
    else
        
        intstatsdummy = [intstatsdummy; intstats(i,:)];
        average_ref_index_dummy = [average_ref_index_dummy; average_ref_index(i)];
        avcangle_dummy = [avcangle_dummy; avcangle(i)];
        cspeed_dummy = [cspeed_dummy; cspeed(i)];
        
        
    end
    
end

intstats = intstatsdummy;
average_ref_index = average_ref_index_dummy;
avcangle = avcangle_dummy;
cspeed = cspeed_dummy;

emptyCells = cellfun('isempty', Detector_Data);
Detector_Data(emptyCells) = [];

newnumint = length(Detector_Data(:,1));


%% Cherenkov angles for each Interaction


% Find deviation from Cherenkov angle for each detector

% Loop each interaction

for i = 1:newnumint
    
    DData = Detector_Data{i,1};
    interaction_point = intstats(i,1:3);
    rc = ring_centre(i,:);
    interactionrcvector = rc - interaction_point; % this is the vector going from the interaction point to the ring centre
    
    %    Loop detectors
    for j = 1:length(DData(:,1))
        
        interactiondvector = DData(j,1:3)- interaction_point; % this is the vector going from the interaction point to the detector
        theta = acos(dot(interactiondvector',interactionrcvector')/(norm(interactionrcvector)*norm(interactiondvector)))-(3.141592653/2); % angle between them (range shifted from default)
        delta_theta = abs(theta) - avcangle(i);
        if (abs(theta) + avcangle(i)) < delta_theta
            delta_theta = abs(theta) + avcangle(i);
        end
        
        DData(j,7) = delta_theta;
        
    end
    
    Detector_Data{i,1} = DData;
end




%% Find E-Field at each detector

Elpm = 2*10^15; % in eV

efieldthresh = 1e-7; %threshold
efieldthresh = sqrt(noise_pow_spec_density * sn_ratio);


% Loop each interaction


for i = 1:newnumint
    
    DData = Detector_Data{i,1};
    DData2 = [];
    interaction_point = intstats(i,1:3);
    attenuation_length = 10*avAttenLength((27-interaction_point(3))*0.1,0.3,attencoeff);
    
    %    Loop detectors
    for j = 1:length(DData(:,1))
        
        Cone_field = e_cherenkov_times_R(i);
        delta_theta = radtodeg(DData(j,7));
        interactiondvector = DData(j,1:3)- interaction_point; % this is the vector going from the interaction point to the detector
        R = norm(interactiondvector)/10;
        E0 = showerE(i);
        dt = 2.7 * (1/nu_scaled) * (Elpm/((0.14*E0)+Elpm))^0.3;
        exponent = - log(2)*(delta_theta/dt)^2;
        field = Cone_field*exp(exponent)*exp(-attenuation_length/R)/R;
        
        if field >= efieldthresh % Only have data for fields above noise threshold
            
            DData2 = [DData2; [DData(j,:) field]];
            
            
        end
        
    end
    
    Detector_Data{i,1} = DData2;
end

% Remove interactions which are undetected or underdetected

emptyCells = cellfun('isempty', Detector_Data(:,1));

Detector_Data(emptyCells,:) = [];
numintdummy = newnumint;

newnumint = length(Detector_Data(:,1));


True_Detector_Data = {};

for i = 1:newnumint
    
    DData = Detector_Data{i,1};
    
    if length(DData(:,1)) >= 5 % Need 5 data points for multilateration
        
        True_Detector_Data{i,1} = DData;
        True_Detector_Data{i,2} = Detector_Data{i,2};
        
        
    end
    
end

Detector_Data = True_Detector_Data;


emptyCells = cellfun('isempty', Detector_Data);

if sum(emptyCells) > 0
    Detector_Data(emptyCells) = [];
    len = length(Detector_Data)/2;
    Detector_Data = reshape(Detector_Data, round(len), 2);
    newnumint = length(Detector_Data(:,1));
    
else
    
    newnumint = 0;
    display('no detectable interactions');
    return
    
end

allowedints = [];

for i = 1:newnumint
    allowedints = [allowedints Detector_Data{i,2}];
end

intstatsdummy = [];
average_ref_index_dummy = [];
avcangle_dummy = [];
cspeed_dummy = [];


for i = 1:length(allowedints)
    
    n = allowedints(i);
    
    % find correct index
    
    index = 0;
    
    for j = 1:numintdummy
        
        a = intstats(j,8);
        
        if a == n
            
            index = j;
            
        end
        
        
    end
    
    intstatsdummy = [intstatsdummy; intstats(index,:)];
    average_ref_index_dummy = [average_ref_index_dummy; average_ref_index(index)];
    avcangle_dummy = [avcangle_dummy; avcangle(index)];
    cspeed_dummy = [cspeed_dummy; cspeed(index)];
    
end

intstats = intstatsdummy;
average_ref_index = average_ref_index_dummy;
avcangle = avcangle_dummy;
cspeed = cspeed_dummy;




% Plot E-Fields at Detectors

figure(7)
hold on
axis([1 100 1 100 1 detector_depth 0 1])

for i = 1:newnumint
    
    DData = Detector_Data{i};
    X = DData(:,1);
    Y = DData(:,2);
    Z = DData(:,3);
    E = DData(:,8);
    S = -log(E);
    S = S./min(S);
    S = S.^3*20;
    scatter3(X,Y,Z,S,'fill')
    
    
end


%% DTOA information

% Time taken

% loop interactions

for k = 1:newnumint
    
    Det_location = Detector_Data{k,1};
    
    
    %    loop detectors
    for j = 1:size(Det_location(:,1))
        
        distance_from_detector = Det_location(j,4);
        time_taken = (distance_from_detector)./(cspeed(k)); % time for radiation to reach each detector
        Det_location(j,5) = time_taken;
        
    end
    
    Detector_Data{k,1} = Det_location;
    
end


% Time differences

% Loop interactions

for k = 1:newnumint
    
    Det_location = Detector_Data{k,1};
    
    % find minimum time
    
    timeMin = min(Det_location(:,5));
    
    % subtract first time to find time difference of signal (TDOA) at each detector in range. Can use these for tri-lateration
    
    Det_location(:,6) = Det_location(:,5) - timeMin;
    
    Detector_Data{k,1} = Det_location;
    
end

%% POLARISATIONS

% E-field is polarised perpendicular to Cherenkov Cone. To find
% polarisation vector find the closest point on the cone to each detector
% and the normalised vector from from this point to the detector is the
% polarisation

for k = 1:newnumint
    
    DData = Detector_Data{k,1};
    
    % Make Cherenkov Cone
    
    cone = makeConeTrunc(intstats(k,:), avcangle(k), detector_depth);
    
    % Loop detectors
    
    polvectorarray = zeros(length(DData(:,1)),3);
    
    for j = 1:length(DData(:,1))
        
        polvector = [];
        distmin = realmax;
        
        dx = DData(j,1);
        dy = DData(j,2);
        dz = DData(j,3);
        
        % Loop through Cone coordinates
        
        if isempty(cone)
            continue
        end
        
        for i = 1:length(cone(1,:))
            
            cx = cone(1,i);
            cy = cone(2,i);
            cz = cone(3,i);
            
            vec = [(dx-cx), (dy-cy), (dz-cz)];
            
            if norm(vec) < distmin % Check if closest point
                polvector = vec/norm(vec);
                distmin = norm(vec);
            end
            
        end
        
        % Put polarisation into detector
        
        polvectorarray(j,:) = polvector;
        
    end
    
    
    
    
    DData = [DData, polvectorarray];
    
    Detector_Data{k,1} = DData;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         EVENT RECONSTRUCTION                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Multilateration from DTOA information

Calculated_positions = [];

%  Comparing sites, loop through each interaction

for k = 1:newnumint
    
    Detector_location = Detector_Data{k,1};
    
    %     Detector_location = sortrows(Detector_location,6);
    
    a = zeros(size(Detector_location(:,1)));
    b = zeros(size(Detector_location(:,1)));
    c = zeros(size(Detector_location(:,1)));
    d = zeros(size(Detector_location(:,1)));
    
    
    for i = 2:size(Detector_location(:,1))
        
        a(i) = (2./cspeed(k)).*((Detector_location(i,1)./(Detector_location(i,6)))-((Detector_location(1,1)./Detector_location(1,6))));
        b(i) = (2./cspeed(k)).*((Detector_location(i,2)./(Detector_location(i,6)))-((Detector_location(1,2)./Detector_location(1,6))));
        c(i) = (2./cspeed(k)).*((Detector_location(i,3)./(Detector_location(i,6)))-((Detector_location(1,3)./Detector_location(1,6))));
        d(i) = (cspeed(k)*Detector_location(i,6))-(cspeed(k)*Detector_location(1,6))-((Detector_location(i,1)^2+Detector_location(i,2)^2 + ...
            Detector_location(i,3)^2)./(cspeed(k)*Detector_location(i,6)))+((Detector_location(1,1)^2+Detector_location(1,2)^2 + ...
            Detector_location(1,3)^2)./(cspeed(k)*Detector_location(1,6)));
        
        
        
    end
    
    aa = a(isfinite(a));
    bb = b(isfinite(b));
    cc = c(isfinite(c));
    dd = d(isfinite(d));
    
    ai = aa';
    bi = bb';
    ci = cc';
    di = dd';
    
    C = [ai; bi; ci]';
    B = -di';
    A = round(C\B);
    
    x = A(1);
    y = A(2);
    
    
    Dsort = sortrows(Detector_location,6);
    
    zmax = detector_depth;
    sumdiff1 = realmax;
    z1 = 0;
    
    
    for z = 1:zmax
        differences = zeros(size(Dsort,1),1);
        for i = 2:size(Dsort,1)
            distance = sqrt((Dsort(i,1)-x)^2 + (Dsort(i,2)-y)^2 + (z-Dsort(i,3))^2);
            distance11 = sqrt((Dsort(1,1)-x)^2 + (Dsort(1,2)-y)^2 + (z-Dsort(1,3))^2);
            differences(i,1) = (distance/cspeed(k)) - (distance11/cspeed(k));
            
        end
        differences = abs(Dsort(:,6)-differences(:,1));
        
        if sumdiff1 > sum(differences)
            z1 = z;
            sumdiff1 = sum(differences);
        end
    end
    
    A(3) = z1;
    A = A'; %gives x,y co-ordinates of interaction
    Calculated_positions = [Calculated_positions; A];
    
    Detector_Data{k,1} = Dsort;
    
end


Calculated_positions %gives output of multilaterated positions

return
%% Velocity Reconstruction

% Loop through velocities, calculate the resulting polarisations at the detectors
% and choose velocity with lowest chi^2 against the measured data.

% Loop interactions

Calculated_velocities = [];

for int = 1:newnumint
    
    h = waitbar(0,'Please wait...');
    DData = Detector_Data{int,1};
    interaction_point = Calculated_positions(int,:);
    
    %    Chi^2 parameters
    right_velocity = [0 0 0];
    
    chi2 = realmax;
    
    % loop velocities
    for i = 1:100 %lower resolution for speed
        waitbar(i/100)
        for j = 1:100 %lower resolution for speed
            for k = 1:detector_depth
                
                vx2 = i-interaction_point(1);
                vy2 = j-interaction_point(2);
                vz2 = k-interaction_point(3);
                interaction_velocity = [vx2 vy2 vz2];
                cone = makeConeTrunc([interaction_point 0 interaction_velocity], avcangle(int), detector_depth);
                if isempty(cone) % Checks if cone doesn't hit surface
                    continue
                end
                % Calculate polarisation
                polvectorarray = zeros(length(DData(:,1)),3);
                
                for jay = 1:length(DData(:,1))
                    
                    polvector = [];
                    distmin = realmax;
                    
                    dx = DData(jay,1);
                    dy = DData(jay,2);
                    dz = DData(jay,3);
                    
                    % Loop through Cone coordinates
                    
                    
                    for ai = 1:length(cone(1,:))
                        
                        cx = cone(1,ai);
                        cy = cone(2,ai);
                        cz = cone(3,ai);
                        
                        vec = [(dx-cx), (dy-cy), (dz-cz)];
                        
                        if norm(vec) < distmin % Check if closest point
                            polvector = vec/norm(vec);
                            distmin = norm(vec);
                        end
                        
                    end
                    
                    % Put polarisation into detector
                    
                    polvectorarray(jay,:) = polvector;
                    
                    
                end
                differences = polvectorarray - DData(:,9:11);
                chi2test = sum(differences.^2); % Add up residuals
                
                if chi2test < chi2
                    chi2 = chi2test;
                    right_velocity = interaction_velocity;
                    
                end
                
                
                
                
                
            end
        end
    end
    
    close(h)
    
    Calculated_velocities = [Calculated_velocities; right_velocity];
    
end


%% Energy Reconstruction

% Loop through energies, calculate resulting E-Fields at relevant detectors
% and choose energy with lowest chi^2 at relevant detectors


% Create dummy detector array

Detector_Data_Dummy = Detector_Data;

% Shorten eV vector

eV2 = [];

for eVindex = 1:length(eV)
    
    enint = eVindex/1000;
    
    if mod(enint,1) == 0
        eV2(enint) = eV(enint);
    end
    
    
end


% Store Energies

Calculated_energies = [];


% Loop interactions
for int = 1:newnumint
    
    
    interaction_point = Calculated_positions(int,:);
    %    Chi^2 parameters
    
    Energy = 0;
    chi2 = realmax;
    
    % Loop velocities
    
    h = waitbar(0,'Please wait...');
    
    
    
    % Find ring centre for this velocity
    
    interaction_velocity = Calculated_velocities(int,:);
    
    lambda = (detector_depth - interaction_point(3))/interaction_velocity(3);
    
    if lambda >= 0 % Upward moving neutrinos
        rx = interaction_point(1) + lambda*interaction_velocity(1);
        ry = interaction_point(2) + lambda*interaction_velocity(2);
        rc = [rx ry detector_depth];
    else % Downward moving neutrinos
        lambda = (1 - interaction_point(3))/interaction_velocity(3);
        rx = interaction_point(1) + lambda*interaction_velocity(1);
        ry = interaction_point(2) + lambda*interaction_velocity(2);
        rc = [rx ry 1];
        
        % Check if ring doesnt hit surface
        vec2check = rc - interaction_point;
        % Dot product with z direction
        dp = dot([0 0 1], vec2check);
        angle = acos(dp)/norm(vec2check)-(3.141592/2);
        if angle > avcangle(int)
            rc = NaN;
        end
    end
    
    if isnan(rc)
        continue
    end
    
    % Find deviations
    
    DData = Detector_Data{int,1};
    DData2 =[];
    
    
    ssddata = sum(size(DData));
    
    if ssddata==0
        
        continue
    end
    
    for col1 = 1:length(DData(:,1))
        
        interactiondvector = DData(col1,1:3)- interaction_point; % this is the vector going from the interaction point to the detector
        theta = acos(dot(interactiondvector',interactionrcvector')/(norm(interactionrcvector)*norm(interactiondvector)))-(3.141592653/2); % angle between them (range shifted from default)
        delta_theta = abs(theta) - avcangle(int);
        if (abs(theta) + avcangle(int)) < delta_theta
            delta_theta = abs(theta) + avcangle(int);
        end
        
        DData(col1,7) = delta_theta;
        
    end
    
    ssddata = sum(size(DData));
    
    if ssddata==0
        
        continue
    end
    
    % Loop Energies
    
    for eVindex = 1:length(eV2)
        E0 = eV2(eVindex);
        
        % Get E-fields
        
        for col2 = 1:size(DData(:,1))
            % Calculate Fields at detectors
            Cone_field = 1.1*10^(-7)*(E0/(1*10^12))*nu_scaled*(1/(1+0.4*nu_scaled^2));
            delta_theta = radtodeg(DData(col2,7));
            interactiondvector = DData(col2,1:3)- interaction_point; % this is the vector going from the interaction point to the detector
            R = norm(interactiondvector)/10;
            dt = 2.7 * (1/nu_scaled) * (Elpm/((0.14*E0)+Elpm))^0.3;
            exponent = - log(2)*(delta_theta/dt)^2;
            DData(col2,8) = Cone_field*exp(exponent)/R;
            
            
            
            if field >= efieldthresh % Only have data for fields above noise threshold
                
                DData2 = [DData2; DData(col2,:)];
                
                
            end
        end
        
        DData = DData2;
        % Now chi^2 with the actual data
        DDatacompare = Detector_Data{int,1};
        
        ssff = sum(size(DData)~=size(DDatacompare)); % Check arrays same size
        if ssff==2
            
            continue
            
        end
        
        deviations = DData(:,8)-DDatacompare(:,8);
        sqdeviations = deviations.^2;
        testchi2 = sum(sqdeviations);
        
        if testchi2 < chi2
            
            chi2 = testchi2;
            
            Energy = E0;
            
        end
        
    end
    
    close(h)
    
    Calculated_energies = [Calculated_energies; E0];
    
end

% Display Laterated Output

detoutput = [Calculated_positions, Calculated_velocities, Calculated_Energies]
