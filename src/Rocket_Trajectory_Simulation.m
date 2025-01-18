%% RPL Launch Script
% Assumptions
% 1. Vertical Takeoff
% 2. Windspeed can be set as constant in horizontal direction
% 3. Air Density is interpolated using natural cubic spline
% 4. Parachute deploys instantaneously
clear;close all;clc;

%% Declare all variables
g=9.81; % m/s

Altitudes=0:10000;
[T, a, P, rho]=atmosisa(Altitudes);

AirDensityFromAltitude = @(x) interp1(Altitudes, rho, x, "linear");

%% Mach Number vs Drag Coefficient Data
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:B314";

% Specify column names and types
opts.VariableNames = ["MachNum", "DragCoef"];
opts.VariableTypes = ["double", "double"];

% Import the data
MachvsDragcoef = readtable("Mach vs Drag coef new.xlsx", opts, "UseExcel", false);

Mach_Num = table2array(MachvsDragcoef(1:300,1));
Drag_coef = table2array(MachvsDragcoef(1:300,2));

% Clear temporary variables
clear opts MachvsDragcoef;

%% Rocket body parameters
DryMass = 54.4218; %kg  (126 lbs)
FuelMass = 20.72917; %kg  (45.7 lbs)
Rocket_Diameter = 0.206756; %m  (8.14 in)
Rocket_Height = 4.64; %m (165 in)
Rocket_Cd = @(v) interp1(Mach_Num.*343, Drag_coef, v, "spline");
Rocket_CrossSectionArea = pi*(Rocket_Diameter/2)^2;

%% Parachute parameters
Drogue_Diameter = 1.524; %m  (5 ft)
Main_Diameter = 4.8768; %m (192 in)
Drogue_Cd = 2.2;
Main_Cd = 2.2;
Main_Deployment_Altitude = 304.8; %m (1000 ft)

Drogue_Area = pi*(Drogue_Diameter/2)^2; %m^2
Main_Area = pi*(Main_Diameter/2)^2; %m^2

%% Engine
m_dot = 1.2579; % kg/s Assume constant mass flow rate
BurnTime= FuelMass/m_dot; %s
mass=@(t) DryMass+FuelMass-m_dot*t; % calculate mass when burning fuel

%% Load Thrust Curve
load("thrust_curve.mat")
Thrust = @(y) interp1(h, T, y, "linear");

%% Windload
avg_Windspeed = 0; %m/s (average 2022 June FAR Site Data: 11.7 mph)
% Windspeed = @(t) (10 + 3*rand(1))*sin(t);
Rocket_WindloadArea = Rocket_Diameter*Rocket_Height;
WindLoad = @(y) 1/2*AirDensityFromAltitude(y)*(avg_Windspeed+randi([-3,3],1,1))^2*Rocket_WindloadArea;
WindLoad_data = (-1 + 2.*rand(10001,1)) .* WindLoad(0:1:10000)';
WindLoad = @(y) interp1(0:1:10000, WindLoad_data, y, "linear");
%% Drag Force
RocketDrag = @(y, vy) 1/2*AirDensityFromAltitude(y).*vy.^2.*Rocket_Cd(vy)*Rocket_CrossSectionArea;
DrogueDrag = @(y, vy) 1/2*AirDensityFromAltitude(y).*vy.^2*Drogue_Cd*Drogue_Area;
MainDrag = @(y, vy) 1/2*AirDensityFromAltitude(y).*vy.^2*Main_Cd*Main_Area;
%% Define System of ODEs to solve
% v = y_dot
% v_dot = 1/m*(F-m_dot*y_dot)
% Solve using ode45

% initial conditions
% q0 = [0;0;0;0];

%% ODE solver
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);

%% Stage 1: Lift off, Engine in operation
disp('Stage 1: Lift off, Engine in operation...');
DE_1 = @(t,q) [q(2);
               % 1/mass(t)*(WindLoad(q(3)) - m_dot*q(2));
               1/mass(t)*WindLoad(q(3));
               q(4);
               % 1/mass(t)*(Thrust(q(3)) - mass(t)*g - RocketDrag(q(3), q(4)) - m_dot*q(4))
               1/mass(t)*(Thrust(q(3)) - mass(t)*g - RocketDrag(q(3), q(4)))];

[t1,stage1] = ode45(DE_1, [0 BurnTime], [0;0;0;0]);

EngineTurnsOff = length(t1);

imp = Thrust(stage1(2:end,3)).*diff(t1);
tot_imp = sum(imp);
Isp = Thrust(0)/(m_dot*g);
fprintf('Total Impulse = %f N*sec \n', tot_imp);
fprintf('Isp at sea level = %f sec \n', Isp);

Force_y = Thrust(stage1(:,3)) - mass(t1)*g - RocketDrag(stage1(:,3), stage1(:,4));
%% Stage 2: Engine turns off, continue toward Apogee
disp('Stage 2: Engine turns off, continuing upward...');

opts2 = odeset('RelTol',1e-5,'AbsTol',1e-7, 'Events', @stop2);

DE_2 = @(t,q) [q(2);
               1/DryMass*WindLoad(q(3));
               q(4);
               1/DryMass*(-DryMass*g - RocketDrag(q(3), q(4)))];

[t2,stage2] = ode45(DE_2, [BurnTime BurnTime+30], [stage1(end,1), stage1(end,2), stage1(end,3), stage1(end,4)], opts2);

Apogee = length(t1)+length(t2);
DrogueDeploy=Apogee;

Force_y = [Force_y; -DryMass*g - RocketDrag(stage2(:,3), stage2(:,4))];
%% Stage 3: Drogue Parachute Opens at Apogee
disp('Stage 3: Drogue parachute deploys at Apogee...');

opts3 = odeset('RelTol',1e-5,'AbsTol',1e-7, 'Events', @stop3);

DE_3 = @(t,q) [q(2);
               1/DryMass*WindLoad(q(3));
               q(4);
               1/DryMass*(RocketDrag(q(3), q(4))+DrogueDrag(q(3), q(4)) - DryMass*g )];

[t3,stage3] = ode45(DE_3, [t2(end) t2(end)+400], [stage2(end,1), stage2(end,2), stage2(end,3), stage2(end,4)], opts3);

MainDeploy = length(t1)+length(t2)+length(t3);

Force_y = [Force_y; RocketDrag(stage3(:,3), stage3(:,4))+DrogueDrag(stage3(:,3), stage3(:,4)) - DryMass*g ];

%% Stage 4: Main parachute deployment + Landing
disp('Stage 4: Main parachute deploys...');

opts4 = odeset('RelTol',1e-5,'AbsTol',1e-7, 'Events', @stop4);

DE_4 = @(t,q) [q(2);
               1/DryMass*WindLoad(q(3));
               q(4);
               1/DryMass*(RocketDrag(q(3), q(4)) + MainDrag(q(3), q(4)) - DryMass*g )];

[t4,stage4] = ode45(DE_4, [t3(end) t3(end)+600], [stage3(end,1), stage3(end,2), stage3(end,3), stage3(end,4)], opts4);

Landing = length(t1)+length(t2)+length(t3)+length(t4);

Force_y = [Force_y; RocketDrag(stage4(:,3), stage4(:,4)) + MainDrag(stage4(:,3), stage4(:,4)) - DryMass*g ];
%% Results

Total_Vel = @(vx,vy) sqrt(vx.^2+vy.^2);

% display some key results
disp('----------------------------------------');
disp('Key Results:');
fprintf('Average wind speed: %.2f m/s\n', avg_Windspeed);
fprintf('Apogee: %f meters\n', stage2(end, 3));
ind = find(isnan(stage4(:,1)));
fprintf('x displacement: %f meters \n', stage4(ind(1)-1, 1));
fprintf('Max Velocity: %f meters/sec\n', Total_Vel(stage2(1,2), stage2(1,4)));
fprintf('Landing velocity: %f meters/sec\n', Total_Vel(stage4(ind(1)-1,2), stage4(ind(1)-1,4)));
fprintf('Total Flight Time: %f sec\n', t4(end));
fprintf('Total Impulse: %f N*sec \n', tot_imp);

%% plot Altitude vs. Time
t = [t1;t2;t3;t4];
y = [stage1(:,3);stage2(:,3);stage3(:,3);stage4(:,3)];
vy = [stage1(:,4);stage2(:,4);stage3(:,4);stage4(:,4)];
x = [stage1(:,1);stage2(:,1);stage3(:,1);stage4(:,1)];
vx = [stage1(:,2);stage2(:,2);stage3(:,2);stage4(:,2)];

figure
hold on;
title('Altitude vs Time')
xlabel('time (sec)')
ylabel('altitude (meters)')
plot(t,y, DisplayName='Altitude vs time');

plot(t(1), y(1), 'r*', DisplayName='Takeoff');
plot(t(EngineTurnsOff), y(EngineTurnsOff), 'r*', DisplayName='EngineTurnsOff');
plot(t(Apogee), y(Apogee), 'r*', DisplayName='Apogee');
plot(t(DrogueDeploy), y(DrogueDeploy), 'r*', DisplayName='Drogue Deployment');
plot(t(MainDeploy), y(MainDeploy), 'r*', DisplayName='Main Deployment');
plot(t(Landing), y(Landing), 'r*', DisplayName='Landing');
legend

%% plot y Velocity vs. Time
figure
hold on;
title('y Velocity vs Time')
xlabel('time (sec)')
ylabel('y velocity (m/s)')
plot(t,vy, DisplayName='Total Velocity vs time');
yline(343, 'g-', DisplayName='Mach 1');

plot(t(1), vy(1), 'r*', DisplayName='Takeoff');
plot(t(EngineTurnsOff), vy(EngineTurnsOff), 'r*', DisplayName='EngineTurnsOff');
plot(t(Apogee), vy(Apogee), 'r*', DisplayName='Apogee');
plot(t(DrogueDeploy), vy(DrogueDeploy), 'r*', DisplayName='Drogue Deployment');
plot(t(MainDeploy), vy(MainDeploy), 'r*', DisplayName='Main Deployment');
plot(t(Landing), vy(Landing), 'r*', DisplayName='Landing');
legend

%% plot Force_y vs. Time
figure
hold on
title('Force on Rocket y vs Time')
xlabel('time (sec)')
ylabel('Force on Rocket y (N)')
plot(t, Force_y, DisplayName='Force on Rocket y vs time');

plot(t(1), Force_y(1), 'r*', DisplayName='Takeoff');
plot(t(EngineTurnsOff), Force_y(EngineTurnsOff), 'r*', DisplayName='EngineTurnsOff');
plot(t(Apogee), Force_y(Apogee), 'r*', DisplayName='Apogee');
plot(t(DrogueDeploy), Force_y(DrogueDeploy), 'r*', DisplayName='Drogue Deployment');
plot(t(MainDeploy+1), Force_y(MainDeploy+1), 'r*', DisplayName='Main Deployment');
plot(t(Landing), Force_y(Landing), 'r*', DisplayName='Landing');
legend

%% plot x displacement from lauching site
figure
hold on;
title('Altitude vs. x-Displacement');
xlabel('x-Displacement (m)');
ylabel('Altitude (m)');
plot(x,y, DisplayName='Altitude vs. Displacement');
plot(x(1), y(1), 'r*', DisplayName='Takeoff');
plot(x(Apogee), y(Apogee), 'r*', DisplayName='EngineTurnsOff');
plot(x(Apogee), y(Apogee), 'r*', DisplayName='Apogee');
plot(x(DrogueDeploy), y(DrogueDeploy), 'r*', DisplayName='Drogue Deployment');
plot(x(MainDeploy), y(MainDeploy), 'r*', DisplayName='Main Deployment');
plot(stage4(end, 1), y(Landing), 'r*', DisplayName='Landing');
legend

%% Define each stop condition
function [value, isterminal, direction] = stop2(~, Y)
    value      = (Y(4) <= 0.01);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end

function [value, isterminal, direction] = stop3(~, Y)
    Main_Deployment_Altitude = 304.8; %m (1000 ft);
    value      = (Y(3) <= Main_Deployment_Altitude);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end

function [value, isterminal, direction] = stop4(~, Y)
    value      = (Y(3) <= 0.01);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end