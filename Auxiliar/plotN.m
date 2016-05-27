function plotN(x)
global g0 Isp
e = [0.001 0.01 0.05 0.1 0.2 0.5];
l = 0.0001:0.0001:1;


figure(2)
for i=1:size(e,2),
U=g0*Isp*log((1+l)./(l+e(i)));
semilogx(l,U)
    hold on
end
U=g0*Isp*log((1+x(1))./(x(1)+x(2)));
semilogx(x(1),U,'b*')
%legend('toggle')
legend('e = 0.001','e = 0.01','e = 0.05','e = 0.1','e = 0.2','e = 0.5')
end
