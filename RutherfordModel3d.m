%RutherfordModel3d.m
%Written by Patrick Hatch 2015/12/07
%University of Western Ontario

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = 1.602e-19; %elementary charge in C
epsilon = 8.854e-12; %permittivity of free space
Z1 = 4; %atomic number for helium (alpha particle)
Z2 = 79; %atomic number for gold (target)
k = 4.87; %energy in MeV
K = k*10^6*e; %conversion of energy into Joules
F = (Z1*Z2*e^2)/(8*pi*epsilon*K); %Important combination

numAlpha = 1000000; %Number of alpha particle
thick = 600*10^-9; %Thickness of foil in m
n = 5.9e28; %number density of atoms (in atoms/m^3)
distAtom = sqrt(1/(n*thick))*10^9; %linear distance between atoms in nm
dFoil = 100; %distance from detector in nm
l = 20*distAtom;  %half length of the foil in nm
thetaRange = 0:0.1:180; %In the following line we set the theoretical curve
radRange = pi*thetaRange/180;
theoScat = (numAlpha*n*thick*(F^2))./(4*((sin(radRange/2)).^4));

alphaMax = atan(l/dFoil); %Max angle a particle can have and still hit the foil
alphaMin = atan(-l/dFoil); %Min angle a particle can have and still hit the foil

%Set up atom positions
m = ceil((2*l)/distAtom);
atomPositions = zeros(1,m);

for jj = 1:m
    atomPositions(jj) = l - distAtom*(jj-1);
end

alphaAngles = zeros(1,numAlpha);

for ii = 1:numAlpha
    %Compute random phi and theta
    phi = 2*alphaMax*rand + alphaMin;
    theta = 2*alphaMax*rand + (pi/2 + alphaMin);
    
    %Compute corrdinates
    x = cos(phi)*sin(theta);
    y = sin(phi)*sin(theta);
    z = cos(theta);
    
    
    %To stretch the vector out such that its x-component is the distance
    %to the foil simply consider that tan(theta) = y/xfoil and rearrange
    %for y
    zfoil = dFoil*(z/x);
    yfoil = dFoil*(y/x);
    xfoil = dFoil;
    
    %Find the closest atom and distance from it
    ymin = min(abs(yfoil - atomPositions));
    zmin = min(abs(zfoil - atomPositions));
    d = sqrt((ymin)^2 + (zmin)^2);
    fromXAxis = acos(x); %Angle from the x-axis
    %Compute impact parameter
    b = d*cos(fromXAxis)*10^-9; %Here we convert the distance into metres for correct units
    
    %Compute scattering angle
    scatterAngle = 2*atan(F/b);
    alphaAngles(1,ii) = scatterAngle;
    
end


%See RutherfordModel.m for reference
alphaAnglesDegrees = 180*alphaAngles/pi;
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

%[n, angleOut] = hist(alphaAnglesDegrees(1,:),0:10:max(alphaAnglesDegrees(1,:)));
%bar(angleOut,n,'barwidth', 2.5 , 'basevalue',10^-1);
%set(gca,'YScale','log')
xlim([0 180])
ylim([1,numAlpha])
hold on
%Plot the theoretical curve
plot(thetaRange,theoScat,'r-');
xlabel('Total Scattering Angle (Degrees)')
ylabel('Number of \alpha Particles')
title('The Scattering Distribution for the Geiger-Marsden Experiment (Rutherford''s Model, Gold,3D)')
legend('Monto Carlo results','Theoretical curve')