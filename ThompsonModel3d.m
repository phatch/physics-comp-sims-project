%ThompsonModel3d.m
%Written by Patrick Hatch 2015/12/07
%University of Western Ontario

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numAlpha = 10000; %Number of alpha particle
dFoil = 100e6; %distance from detector in nm
l = 10e6;  %half length of the foil in nm
thick = 600; %thickness of the foil in nm (can be a vector for 
%investigating multiple thicknesses)
distAtom = 0.26; %linear distance between atoms in nm

maxangle = zeros(1,length(thick));

me = 0.511; %MeV/c^2
mAlpha = 3727; %MeV/c^2

mRatio = me/mAlpha;

scatterMax = atan(2*mRatio); %Max deflection
scatterMin = atan(-2*mRatio); %Min deflection

alphaMax = atan(l/dFoil); %Max angle a particle can have and still hit the foil
alphaMin = atan(-l/dFoil); %Min angle a particle can have and still hit the foil


thickcount = 1;

for kk = thick %This loop is only necessary for investigating multiple
    %thickness
    diffAngles = zeros(1,numAlpha);

for ii = 1:numAlpha
    %Random angles
    phi = 2*alphaMax*rand + alphaMin;
    theta = 2*(alphaMax)*rand + (pi/2 + alphaMin);
    
    %Resulting coordinate
    x = cos(phi)*sin(theta);
    y = sin(phi)*sin(theta);
    z = cos(theta);
    rInitial = [x,y,z];
    
    %To stretch the vector out such that its x-component is the distance
    %to the foil simply consider that tan(theta) = y/xfoil and rearrange
    %for y (same for z)
    zfoil = dFoil*(z/x);
    yfoil = dFoil*(y/x);
    xfoil = dFoil;

    
    while 1 %Loop until an if statement breaks it
        %Begin the random walk by computing a random angle between
        %a limited range 
        scattertheta = 2*scatterMax*rand + scatterMin;
        scatterphi = 2*scatterMax*rand + scatterMin;
        oldtheta = theta;
        oldphi = phi;
        theta = theta + scattertheta;
        phi = phi + scatterphi;
        
        rAdd = distAtom*([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)]);
        
        %Uncomment following lines to restrict angles withing a cone of
        %0.016 degrees
        
        %newDisp = [xfoil,yfoil,zfoil] + rAdd;
        
        %proj = acos(dot([xfoil,yfoil,zfoil],newDisp)/(norm([xfoil,yfoil,zfoil])*norm(newDisp)));
        
        %if abs(proj) > scatterMax %If the projection of our displacement vector
            %along the original vector is greater than the maximum amount
            %of scattering allowed by the Thompson model, then we try again
            %theta = oldtheta;
            %phi = oldphi;
            %continue
        %end
            
        %Then compute the resulting direction vector
        %and add it to the previous position
        xfoil = xfoil + rAdd(1);
        yfoil = yfoil + rAdd(2);
        zfoil = zfoil + rAdd(3);
        
        %if the particle exits the bounds of the foil, stop
        if xfoil > (dFoil + kk) || xfoil < dFoil || abs(yfoil) > l || abs(zfoil) > l
            %Compute final scattering angle
            w = dot(rInitial,rAdd);
            a = w/distAtom;
            projFinal = acos(a);
            diffAngles(ii) = projFinal;
            
            break
        end
            
    end
    
end
%Record the max thickness
maxangle(thickcount) = max(diffAngles);
thickcount = thickcount + 1;
end


%Uncomment the lines below when investigating multiple thickness
%maxangleDegrees = 180*maxangle/pi;
%p = polyfit(thick,maxangleDegrees.^2,1);
%thickRange = 100:20000;
%fit = polyval(p,thickRange);
%plot(thickRange,fit)
%hold on
%plot(thick,maxangleDegrees.^2,'ro')
%xlabel('Gold-foil thinkness (nm)')
%ylabel('Max scattering angle (degrees)')
%title('Foil thickness versus max scattering angle')

%Plot histogram, comment out when investigating multiple thickness
%i.e. Figure 5
diffAnglesDegrees = 180*diffAngles/pi;


hist(diffAnglesDegrees,40);
xlabel('Total Scattering Angle (Degrees)')
ylabel('Number of \alpha Particles')
title('The Scattering Distribution for the Geiger-Marsden Experiment (Thompson''s model,Gold,3D)')
