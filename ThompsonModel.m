%ThompsonModel.m
%Written by Patrick Hatch 2015/12/06
%University of Western Ontario

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numAlpha = 10000; %Number of alpha particle
dFoil = 10e6; %distance from detector in nm
l = 1e6;  %half length of the foil in nm
thick = 600; %thickness of the foil in nm (can be a vector for 
%investigating multiple thicknesses 
distAtom = 0.26; %linear distance between atoms in nm

maxangle = zeros(1,length(thick));

me = 0.511; %MeV/c^2
mAlpha = 3727; %MeV/c^2

mRatio = me/mAlpha;

scatterMax = atan(2*mRatio);
scatterMin = atan(-2*mRatio);

thetaMax = atan(l/dFoil); %Max angle a particle can have and still hit the foil
thetaMin = atan(-l/dFoil); %Min angle a particle can have and still hit the foil

projAngles = zeros(1,numAlpha);
thickcount = 1;

for kk = thick %This loop is only necessary for investigating multiple
    %thickness
    alphaAngles = zeros(2,numAlpha);
    diffAngles = zeros(1,numAlpha);
for ii = 1:numAlpha
    %Random angle
    theta = 2*thetaMax*rand + thetaMin;
    alphaAngles(1,ii) = theta;
    
    %Resulting coordinate
    x = cos(theta);
    y = sin(theta);
    rInitial = [x,y];
    
    %To stretch the vector out such that its x-component is the distance
    %to the foil simply consider that tan(theta) = y/xfoil and rearrange
    %for y
    yfoil = dFoil*(y/x);
    xfoil = dFoil;

    
    while 1 %Loop until an if statement breaks it
        
        %Begin the random walk by computing a random angle between
        %a limited range 
        scatterAngle = 2*scatterMax*rand + scatterMin;
        theta = theta + scatterAngle;
        
        %Then compute the resulting direction vector
        %and add it to the previous position
        rAdd = distAtom*[cos(theta),sin(theta)];
        xfoil = xfoil + rAdd(1);
        yfoil = yfoil + rAdd(2);
        
        %if the particle exits the bounds of the foil, stop
        if xfoil > (dFoil + kk) || xfoil < dFoil || abs(yfoil) > l
            alphaAngles(2,ii) = theta;
            %Compute final scattering angle
            diffAngles(ii) = abs(alphaAngles(1,ii) - alphaAngles(2,ii));
            
            break
        end
            
    end
    
    
end
%Record the max thickness
maxangle(thickcount) = max(diffAngles);
thickcount = thickcount + 1;
end

maxangleDegrees = 180*maxangle/pi;

%Uncomment the lines below when investigating multiple thickness
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
title('The Scattering Distribution for the Geiger-Marsden Experiment (Thompson'' Model, Gold,2D)')
    
    
    
    
    
        
        
        
    
    
    
    
    
    