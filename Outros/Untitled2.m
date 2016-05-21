close all
clear all
clc

Isp = 282;
m0 = 541300;
mp = 507500;
mp1 = 398887;
mp2 = 108613;
m = 196583.0189;
mdot = 273.5849057;
mdot1 = 273.585048;
mdot2 = 273.5843829;
g(1) = 9.81;
v(1) = 0;
y(1) = 0;
m(1) = m0;
t = 0;
for cont = 2:1:140,
t = t + 1;
m(cont) = m0 - mdot*t*9;
v(cont) = v(cont-1) + 3.6*g(cont-1)*(-Isp*log(m(cont)/m0)-t);
y(cont) = y(cont-1) + g(cont-1)*(-t*Isp*((log(m0/m(cont)))/((m0/m(cont))-1))+t*Isp-0.5*t*t);
%g(t+1) = 9.81/((6400+y(t))/6400)^2
g(cont) = 9.81/((6400+y(cont))/6400)^2;
end
y = y(:)/1000;
t = 0:1:cont-1;
figure(1)
plot(t,y,'r*')
xlabel('')
ylabel('Pos [km]')
figure(2)
plot(t,v,'r*')
xlabel('tempo (s)')
ylabel('Vel. [km/s]')
figure(3)
plot(t,g,'r*')
xlabel('tempo (s)')
ylabel('Gravidade [m/ss]')
