clear; close all; clc
m=200;
u=8000;
q=1;
g=32.17;
b=100;
 
count=1;
h(1)=0;
for t=0:3:300
   
    if t<=b
        v(count)=u*log(m/(m-q*t))-g*t;
    else  
        v(count)=u*log(m/(m-q*b))-g*t;
    end
   
 
    h(count+1)=h(count)+v(count)*3;
   
    if h(count+1)<0;
        h(count+1)=0;
    end
   
    tt(count)=t;
    count=count+1;
   
end
 
figure(1)
subplot(121)
plot(tt,v)
xlabel('Time (sec) ');
ylabel('Velocity (ft/sec) ');
title('Rocket Velocity Curve ');
grid
axis tight
 
subplot(122)
N=length(tt);
% h array is one longer than tt.
plot(tt,h(1:N))
xlabel('Time (sec) ');
ylabel('Height (ft) ');
title('Rocket Height Curve ');
grid
axis tight
 
     % display the max height and when that happens in command window
[maxh,nmax]=max(h);
fprintf('Rocket reaches %.1f feet after %.1f seconds \n',maxh,tt(nmax));
 