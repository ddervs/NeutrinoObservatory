%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  	INTERACTION LOCATIONS                  		 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% This script generates a 10x10x27 (or 28) array with each element containing
% a corresponding probability that a given neutrino interaction will take
% place in that position. The code takes into account the possible paths a
% neutrino can take through space, the earth and the ice to get to a given
% position.


wa = waitbar(0,'Progress');
nrow = 10;
ncol = 10;
ndep = 27; % 28 for Antarctica




normalised_probabilities = zeros(nrow,ncol,ndep);
unnormalised_probabilities = zeros(nrow,ncol,ndep);


ecentre = -6371+ndep;
eradius = 6371;
c = [(nrow/2) (ncol/2) ecentre];
[xs,ys,zs] = sphere(25);
xs = reshape(xs.*eradius+5,676,1);
ys = reshape(ys.*eradius+5,676,1);
zs = reshape(zs.*eradius+ecentre,676,1);
XYZ = [xs,ys,zs];
XYZ = sortrows(XYZ,3);
XYZ = XYZ(26:650,:); % Remove top and duplicates

scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3))
interaction_length_earth = 1.66*10^-27/(10^-33*5515);


%% Density profile

den = 25.75; % Antarctica: Use den = 300
gam = 0.022; % 		       gam = 0.04	
AA = 1.1;    %		       AA = 325	
bet = 75;    %		       bet = 1	 	 	

denscoeff = [den AA gam bet];

%% Calculation

for z0 = 1:ndep;
    
    for x0 = 1:nrow;
        for y0 = 1:ncol;
            
            k = 1; % top face
            for i= 1:ncol
                for j= 1:nrow
                    edistnprob = 0;
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2); % modulus of vector from face to point
                    
                    if lmod ~= 0
                        l1 = (x0-i)/lmod;
                        l2 = (y0-j)/lmod;
                        l3 = (z0-k)/lmod;
         
                        o = [i j k]; % point on face position vector
                        l = [l1 l2 l3]; % unit vector from point to face
                        
                        
                        l_dot_o_minus_c = dot(l, (o-c));
                        o_minus_c_sqrd = dot((o-c),(o-c));
                          
                        discr = sqrt((l_dot_o_minus_c)^2-o_minus_c_sqrd+eradius^2); % check for real solutions
                        
                        if imag(discr) == 0
                            edistn = max(abs([-l_dot_o_minus_c + discr; -l_dot_o_minus_c - discr]));
                            edistnprob = exp(-edistn/interaction_length_earth);
                            
                            
                        end
                        
                    end
                    
                    if z0 ~= k % Dont use same oint as point and face
                        
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff)); % ANTARCTICA: use avRefIndex() as slightly different code (no additive factor)
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + ((1-exp(-lmod*10/interaction_length_ice))*edistnprob);
                    
                    end
                    
                end
            end
            
            k = ndep;
            for i= 1:ncol
                for j= 1:nrow
                    
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2);
                    
                    if z0 ~= k
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff));
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + (1-exp(-lmod*10/interaction_length_ice));
                    end
                end
            end
            
            i = 1;
            
            
            for k= 1:ndep
                for j= 1:nrow
                    edistnprob = 1;
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2);
                    if lmod ~= 0
                        
                        l1 = (x0-i)/lmod;
                        l2 = (y0-j)/lmod;
                        l3 = (z0-k)/lmod;
                        o = [i j k];
                        l = [l1 l2 l3];
                        
                        
                        
                        l_dot_o_minus_c = dot(l, (o-c));
                        o_minus_c_sqrd = dot((o-c),(o-c));
                        
                        discr = sqrt((l_dot_o_minus_c)^2-o_minus_c_sqrd+eradius^2);
                        
                        if imag(discr) == 0
                            edistn = min(abs([-l_dot_o_minus_c + discr; -l_dot_o_minus_c - discr]));
                            edistnprob = exp(-edistn/interaction_length_earth);
                        end
                        
                    end
                    if z0 ~= k
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff)); % ANTARCTICA: use avRefIndex() as slightly different code (no additive factor)
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + ((1-exp(-lmod*10/interaction_length_ice))*edistnprob);
                    end
                end
            end
            
            i=ncol;
            
            for k= 1:ndep
                for j= 1:nrow
                    edistnprob = 1;
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2);
                    
                    if lmod ~= 0
                        l1 = (x0-i)/lmod;
                        l2 = (y0-j)/lmod;
                        l3 = (z0-k)/lmod;
                        o = [i j k];
                        l = [l1 l2 l3];
                        
                        
                        l_dot_o_minus_c = dot(l, (o-c));
                        o_minus_c_sqrd = dot((o-c),(o-c));
                        
                        discr = sqrt((l_dot_o_minus_c)^2-o_minus_c_sqrd+eradius^2);
                       
                        if imag(discr) == 0
                            edistn = min(abs([-l_dot_o_minus_c + discr; -l_dot_o_minus_c - discr]));
                            edistnprob = exp(-edistn/interaction_length_earth);
                        end
                        
                    end
                    if z0 ~= k
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff)); % ANTARCTICA: use avRefIndex() as slightly different code (no additive factor)
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + ((1-exp(-lmod*10/interaction_length_ice))*edistnprob);
                    end
                end
            end
            
            j = 1;
            
            for k= 1:ndep
                for i= 1:ncol
                    edistnprob = 1;
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2);
                    if lmod ~= 0
                        l1 = (x0-i)/lmod;
                        l2 = (y0-j)/lmod;
                        l3 = (z0-k)/lmod;
                        o = [i j k];
                        l = [l1 l2 l3];
                        
                        
                        l_dot_o_minus_c = dot(l, (o-c));
                        o_minus_c_sqrd = dot((o-c),(o-c));
                        
                        discr = sqrt((l_dot_o_minus_c)^2-o_minus_c_sqrd+eradius^2);
                        if imag(discr) == 0
                            edistn = min(abs([-l_dot_o_minus_c + discr; -l_dot_o_minus_c - discr]));
                            edistnprob = exp(-edistn/interaction_length_earth);
                        end
                        
                    end
                    
                    if z0 ~= k
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff)); % ANTARCTICA: use avRefIndex() as slightly different code (no additive factor)
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + ((1-exp(-lmod*10/interaction_length_ice))*edistnprob);
                    end
                end
            end
            
            j = nrow;
            
            for k= 1:ndep
                for i= 1:ncol
                    
                    lmod = sqrt((x0-i)^2+(y0-j)^2+(z0-k)^2);
                    edistnprob = 1;
                    
                    if lmod ~= 0
                        l1 = (x0-i)/lmod;
                        l2 = (y0-j)/lmod;
                        l3 = (z0-k)/lmod;
                        o = [i j k];
                        l = [l1 l2 l3];
                        
                        l_dot_o_minus_c = dot(l, (o-c));
                        o_minus_c_sqrd = dot((o-c),(o-c));
                        
                        discr = sqrt((l_dot_o_minus_c)^2-o_minus_c_sqrd+eradius^2);
                        
                        if imag(discr) == 0
                            edistn = min(abs([-l_dot_o_minus_c + discr; -l_dot_o_minus_c - discr]));
                            edistnprob = exp(-edistn/interaction_length_earth);
                            
                        end
                        
                    end
                    
                    if z0 ~= k
                        interaction_length_ice = 1.66*10^-27/(10^-33*avDensity((100*z0),(100*k),denscoeff)); % ANTARCTICA: use avRefIndex() as slightly different code (no additive factor)
                        
                        normalised_probabilities(x0,y0,z0) =  normalised_probabilities(x0,y0,z0) + ((1-exp(-lmod*10/interaction_length_ice))*edistnprob);
                    end
                    
                end
            end
            
        end
        
        
    end
    
    
    waitbar(z0/ndep,wa)
end

delete(wa);

unnormalised_probabilities = normalised_probabilities;
normalised_probabilities = normalised_probabilities./sum(sum(sum(normalised_probabilities)));


%% PLOT

sperl = normalised_probabilities;
figure(1)
q = slice(sperl, 1:size(sperl,1), 1:size(sperl,2), 1:size(sperl,3));
alpha('color')
set(q,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
xlabel('width of detector [x100m]')
ylabel('length of detector [x100m]')
zlabel('depth of detector [x100m]')
title('probability of an interaction at any location')
