clc; clear;
%% Finding P_SJ 
Ae = 0.0157770468; % Exit Area
A_star = 0.0029186351; % Critical Area
P0 = 25*101325; % Pa
Area_Ratio = Ae/A_star;
[T, a, P, rho] = atmosisa(0:10000);

Pa = P; %Pa
%%  For P_SJ conditon
syms M
y = 1.2 ;%specific heat ratio
R = 283; %Gas Constant
T0 = 3225; %Kelvin
Area_eq = 1./M*(2/(y+1)*(1+(y-1)/2*M.^2))^((y+1)/(2*(y-1))); % Area equation
Me =double(vpasolve(Area_eq == Area_Ratio,M, [0 inf])); % 
Me = Me(Me>1); % Exit Mach Number for supersonic Jet condition
temp = 1+(y-1)/2.*Me.^2; %temporary var
P_SJ = 1./(temp.^(y/(y-1))) * P0; % Supersonic Jet Pressure
SJ_alt = find(abs(Pa-P_SJ)./P_SJ*100 < 0.01, 1); %altitude at which SuperSonic Jet condition is met
fprintf('P_SJ = %g Pa occurs at an approximate altitude of %g meters\n',P_SJ,SJ_alt);

mdot = 1.2579; % mast flow rate (kg/s)
Temp = temp^(-1) *T0; % Exit Temperature (K)
a = sqrt(y*R*Temp); %speed of sound at exit (m/s)
Ue = Me*a; % Exit Velocity (m/s)
%% Save Thrust Curve to mat file
T = mdot.*Ue + (P_SJ-Pa).*Ae;
h = 0:10000; %height (m)
save("thrust_curve.mat", 'h','T');

%% Thrust Curves
figure
plot(Pa/101325,T,'LineWidth',1)
hold on;
plot(Pa(end)/101325,mdot.*Ue + (P_SJ-Pa(end)).*Ae,'*')
legend('','Maximum Thrust (N)')
title('Thrust vs Atmospheric Pressure')
xlabel('Atmospheric Pressure (Pa)')
ylabel('Thrust (N)')
ylim([3200,4600])
hold off;

figure
h = 0:10000; %height (m)
plot(h,T,'linewidth',1);
hold on;
plot(h(end),mdot.*Ue + (P_SJ-Pa(end)).*Ae,'*')
title('Thrust vs Altitude')
legend('','Maximum Thrust (N)')
xlabel('Altitude (m)')
ylabel('Thrust (N)')
ylim([0,6000])
xlim([0,12000])