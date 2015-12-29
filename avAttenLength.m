function [ atten_length ] = avAttenLength( start_depth, end_depth, polycoeff )
 
% avAttenLength finds the average attenuation length along a given depth
% rang

%% ALL UNITS IN km

latb = 0;
lata = 0;
for i = 1:6

    latb = latb + (polycoeff(i)*end_depth^i)/i;
    lata = lata + (polycoeff(i)*start_depth^i)/i;

end

atten_length = (latb - lata)/(10*(end_depth-start_depth));

end

