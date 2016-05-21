clc
clear all

%% liquid
IspSL = 311; %s
IspVAC = 338; %s
tb = 253; %s

%solid
Isp = 330; %s
tbS = 94; %s

mp = 2889; %kg
mdot = 900; %kg/s


m0 = 250000;
m = m0;
g0 = 9.81;
y = 0;
dt = 1;
v = 0
for t = 0:dt:100;
y = [y (g0*IspSL*t*(1 - log(m0/m)/(((m0/m)-1))) - (1/2)*g0*t*t )];
g0*IspSL*t*(1 - log(m0/m)/((m0/m-1)));
v = [v (-g0*IspSL*log(m/m0) - g0*t)];
m = m -mdot*dt;
%(1/2)*g0*t*t
end

t = 0:dt:100;
t = [0 t];
figure(1)
plot(t,y,'*');
figure(2)
plot(t,v,'*');

