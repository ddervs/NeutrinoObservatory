
clear all
clc
numofobservedints = [];
numofintstakenplace = [];
simulits = 10000;
rootN = 1/sqrt(simulits);


for asdasf = 1:simulits


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
        detector_depth = 27; % times 100m

        %% Density profile

        den = 25.75; % Coefficients of density profile
        gam = 0.022;
        AA = 1.1;
        bet = 75;

        denscoeff = [den AA gam bet];

        rho = avDensity(0,100*detector_depth,denscoeff); %gcm-3
        l = u/(xsec*(rho/1000)); %interaction length in cm


        %% Probability of interaction location
        sperl = new_probabilities;
        S = sperl;

%         figure(1) % Plot location likelihoods
%         q = slice(sperl, 1:size(sperl,1), 1:size(sperl,2), 1:size(sperl,3));
%         alpha('color')
%         set(q,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
%         xlabel('width of detector [x100m]')
%         ylabel('length of detector [x100m]')
%         zlabel('depth of detector [x100m]')
%         title('probability of an interaction at any location')
% 
% 
%         hhbar = waitbar(0,'progress');
       % for jbb = 1:1; % Loop needed to show waitbar 

        %% Variables

        % ice density g/cm3 varying with temp
        Rho = rho/1000 + (5.4E-3.*rand([1,n1])); 

        % ice depth (cm) randomly +/-50cm
        d = 270000+((50-(100.*rand([1,n1]))));

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

%         display('Interactions in the detector per year')
        Totalinteractions = mean(Ninttot);

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

        % %% Plot of Neutrino Energy vs Flux
        
%         figure(3)
        
        bf = polyfit(log(eV), log(Tflux), 1);
        % compute fit in linear domain
        Tflux_hat = exp(bf(1) * log(eV) + bf(2));
        % make log log plot
%         loglog(eV, Tflux_hat,'.', eV, Tflux);
%         label = ['log(Tflux) = ' num2str(bf(1)) 'log(eV) + ' num2str(bf(2))];
%         legend('data', label);
        % hold off;
        
        [best_f] = 10.^(bf(1,1)*log(sort(Tflux))+bf(1,2));%will give weightings to energy of interacting neutrino
        
        
        % 
        for tet = 1:(round(Totalinteractions))% # possible numbers
            ce = sort(eV);
            w3 = sort(best_f,'descend'); %
            N = 1;              %# how many numbers to generate
        
            ei(tet) = ce( sum( bsxfun(@ge, rand(N,1), cumsum(w3./sum(w3))), 2) + 1 ); % This vector is the energies of the interactions
        end
        % 



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



        % %% Plot interaction locations and paths through ice to get there
        % 
        % figure(4)
        % xlabel('width of detector [x100m]')
        % ylabel('length of detector [x100m]')
        % zlabel('depth of detector [x100m]')
        % title('Location and path of Interactions + radiation')
        % for ppo = 1:round(Totalinteractions)
        % grid on
        % hold on
        % stats_plotable = [interaction_stats(ppo,1:3);interaction_stats(ppo,5:7)]';
        % plot3([stats_plotable(2,1),stats_plotable(2,2)],[stats_plotable(1,1),stats_plotable(1,2)],[stats_plotable(3,1),stats_plotable(3,2)],'--')
        % 
        %     plot3(stats_plotable(2,1),stats_plotable(1,1),stats_plotable(3,1),'o')
        % 
        % end
        % hold off
        % axis([0 100 0 100 0 27])
        % daspect([1 1 1])
        % waitbar(jbb/10)
        % end
        % delete(hhbar)
        % 
        % interaction_stats(:,5:7) = interaction_stats(:,5:7)-position_of_interactions;
        % interaction_stats
        % display('        ------------------------------------------------------------------------------')
        % display('           x           y           z     v(m/s)            i           j           k ')

        intstats = []; % Matrix to be filled with interaction information only for detected interactions





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                               DETECTION                               %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % This section of the code proceeds by taking the Monte-Carlo generated
        % interactions, and finding which of the interactions would be detected. 
        % The code works via a process of elimination, eventually returning a cell 
        % array, Detector_Data. Each element corresponds to one of the DETECTED 
        % interactions and is an array of doubles. Each row of the array 
        % corresponds to one detector which has received a signal from the neutrino
        % interaction. Along each row, (in order) is: the x,y,z co-ordinates of the
        % detector, the distance from the interaction, the time the signal takes to
        % get to the detector, the time difference from the detector's lowest time 
        % recorded, the angular deviation from the Cherenkov cone and the Magnitude
        % of the E-Field observed at the detector.

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
%         figure(5)
%         hold on
%         daspect ([1 1 1])
%         axis([0 100 0 100 0 detector_depth])
        [x,y,z] = sphere;      %# Makes a 21-by-21 point sphere
        x = x(1:11,:);       %# Keep top 11 x points
        y = y(1:11,:);       %# Keep top 11 y points
        z = z(1:11,:);       %# Keep top 11 z points

        % %% Plot Detection Spheres
        % 
        Detector_location = New_det_loc;
        % for j = 1:size(Detector_location)
        % dsp = mesh((detection_radius.*x)+(Detector_location(j,1)),(detection_radius.*y)+(Detector_location(j,2)),(detection_radius.*z)+(Detector_location(j,3)));  %# Plot the surface
        %     
        % end
        % 
        % for k = 1:numintround
        %         scatter3(position_of_interactions(k,1),position_of_interactions(k,2),position_of_interactions(k,3),'o','fill','k','LineWidth',6)
        % end



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

        %     ddatasize = length(Det_location(:,1));
        %     Det_location = [Det_location, zeros(ddatasize,4)];    
        %                                                        % Empty slots for DTOA information, added here
        %     Det_location(:,9) = k;                             % due to way code was originally structured



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
%             display('no detectable interaction');
            newnumint = numintround;
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

        refind0 = 1.3;
        A = 0.225;
        gamma = 0.037;
        beta = 1.1;

        ncoeff = [refind0 0.225 0.037 1.1];
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

        x6 = -4E-15;
        x5 = 3E-12;
        x4 = 6E-10;
        x3 = 1E-06;
        x2 = 0.0004;
        x1 = 0.0233;
        x0 = - 1.4911;
        attencoeff = [x0 x1 x2 x3 x4 x5 x6];



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

            
%             display('no detectable interactions'); 

        end



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

