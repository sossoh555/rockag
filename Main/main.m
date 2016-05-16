clc
clear all
close all
global m0 g0 T A Cd Re hgr_turn tf md %rh0 H0
% Launch Site: Guiana Space Center
Alt = 1;              %[m] Alt above sea level

% VEGA Rocket
m_stage_gross = [95796, 25751,10948];% 1st, 2nd,3d
% First stage(Solid Fuel)   
m_prop = 507500;      % *GA* {300,000 - 700,000}[kg] Propellant mass
Isp    = 282 ;        % *CONSTANT*[s]  Specific impulse
d      = 3;           % *GA* {2 - 6} [m]  Diameter 
g0     = 9.81;        % [m/s^2] Constant at its sea-level value
m0  = 541300;         % *GA* { [kg] Initial mass
A   = pi*d^2/4;       % [m^2]Frontal area
Cd  = 0.5 ;           % Drag coefficient,assumed to have the constant value
%rh0 = 1.225;          % [kg/m^3]
%H0 = 7500;            % [m] Density scale height
Re = 6378e3;          % [m] Earth's radius
hgr_turn = 200;       % [m] Rocket starts the gravity turn when h = hgr_turn
tburn = 306;          % [s] Fuell burn time, first stage
md = (m_prop)/tburn;  % [kg/s]Propellant mass flow rate
T   = md*(Isp*g0);    % [N] Thrust (mean)
mf = m0 - m_prop;     % [kg] Final mass of the rocket(first stage is empty)
t0 = 0;               % Rocket launch time
tf = t0 + tburn;      % The time when propellant is completely burned
%and the thrust goes to zero
t_range     = [t0,tf];  % Integration interval

% Launch initial conditions:
gamma0 = 89.5/180*pi;       % Initial flight path angle
v0 = 0;   % Velocity (m/s)  % Earth's Rotation considered in eq of motion.
x0 = 0;   % Downrange distance [km]
h0 = Alt; % Launch site altitude [km]
vD0 = 0;  % Loss due to drag (Velocity)[m/s]
vG0 = 0;  % Loss due to gravity (Velocity)[m/s]
state0   = [v0, gamma0, x0, h0, vD0, vG0];
% Solve initial value problem for ordinary differential equations
[t,state] = ode45(@RocketDynEq,t_range,state0) ;
v     = state(:,1)/1000;      % Velocity [km/s]
gamma = state(:,2)*180/pi;    % Flight path angle  [deg]
x     = state(:,3)/1000;      % Downrange distance [km]
h     = state(:,4)/1000;      % Altitude[km]
vD    = -state(:,5)/1000;     % Loss due to drag (Velocity)[m/s]
vG    = -state(:,6)/1000;     % Loss due to gravity (Velocity)[m/s]
plot(t,h,'.b');
hold on;
grid on;
plot(t,h,'.b');
title('Rocket Dynamics');
xlabel('time [s]');
ylabel('Altitude [km]');


% VEGA Rocket: First Stage P80
fprintf('\n VEGA Rocket: First Stage P80\n')
fprintf('\n Propellant mass           = %4.2f [kg]',m_prop)
fprintf('\n Gross mass                = %4.2f [kg]',m_stage_gross(1))
fprintf('\n Isp                       = %4.2f [s]',Isp)
fprintf('\n Thrust(mean)              = %4.2f [kN]',T/1000)
fprintf('\n Initial flight path angle = %4.2f [deg]',gamma0*180/pi)
fprintf('\n Final speed               = %4.2f [km/s]',v(end))
fprintf('\n Final flight path angle   = %4.2f [deg]',gamma(end))
fprintf('\n Maximum Altitude          = %4.2f [km]',max(h(:)))
fprintf('\n Downrange distance        = %4.2f [km]',x(end))
fprintf('\n Drag loss                 = %4.2f [km/s]',vD(end))
fprintf('\n Gravity loss              = %4.2f [km/s]',vG(end))
fprintf('\n Mass Flow                 = %4.2f [kg/s]',md)
fprintf('\n');