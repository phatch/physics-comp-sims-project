%RutherfordModel.m
%Written by Patrick Hatch 2015/12/06
%University of Western Ontario

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = 1.602e-19; %elementary charge in C
epsilon = 8.854e-12; %permittivity of free space
Z1 = 4; %atomic number for helium (alpha particle)
Z2 = 79; %atomic number for gold (target)
k = 6.5; %energy in MeV
K = k*10^6*e; %conversion of energy into Joules
F = (Z1*Z2*e^2)/(8*pi*epsilon*K); %Important combination

numAlpha = 1000000; %Number of alpha particle
thick = 600*10^-9; %Thickness of foil in m
n = 5.9e28; %number density of atoms (in atoms/m^3)
distAtom = sqrt(1/(n*thick))*10^9; %linear distance between atoms in nm
dFoil = 100; %distance from detector in nm
l = 20*distAtom;  %half length of the foil in nm


thetaMax = atan(l/dFoil); %Max angle a particle can have and still hit the foil
thetaMin = atan(-l/dFoil); %Min angle a particle can have and still hit the foil

%Here we equal space atoms throughout the foil
m = ceil((2*l)/distAtom);
atomPositions = zeros(1,m);

for jj = 1:m
    atomPositions(jj) = l - distAtom*(jj-1);
end



for ii = 1:numAlpha %For each alpha particle
    %Assign a random angle
    theta = 2*thetaMax*rand + thetaMin;
    
    %Compute the x and y coordinates
    x = cos(theta);
    y = sin(theta);
    
    %To stretch the vector out such that its x-component is the distance
    %to the foil simply consider that tan(theta) = y/xfoil and rearrange
    %for y
    yfoil = dFoil*(y/x);
    xfoil = dFoil;
    
    %Then find the closest atom
    d = min(abs(yfoil - atomPositions));
    %Then compute the impact parameter
    b = d*cos(theta)*10^-9; %Here we convert the distance into metres for correct units
    %And then find the scattering angle
    scatterAngle = 2*atan(F/b);
    alphaAngles(1,ii) = scatterAngle;
end


alphaAnglesDegrees = 180*alphaAngles/pi;

%This solution to creating a log histogram is courtesy of 
%http://www.mathworks.com/matlabcentral/answers/91621-
%why-does-my-histogram-become-incorrect-when-i-change-the-y-axis-scaling-to-log
hist(alphaAnglesDegrees,15)
ph = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph);
for i = 1:N_patches
      % Get patch vertices
      vn = get(ph(i),'Vertices');
      % Adjust y location
      vn(:,2) = vn(:,2) + 1;
      % Reset data
      set(ph(i),'Vertices',vn)
end
% Change scale
set(gca,'yscale','log')

xlim([0 180])
ylim([1,numAlpha])
xlabel('Total Scattering Angle (Degrees)')
ylabel('Number of \alpha Particles')
title('The Scattering Distribution for the Geiger-Marsden Experiment (Rutherford''s Model, Gold,2D)')




    
    
    
