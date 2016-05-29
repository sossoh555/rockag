function plotN(x)
global g0 Isp cut test type out finalName PATH
e = [0.001 0.01 0.05 0.1 0.2 0.5];
l = 0.0001:0.0001:1;

if strcmp(type,'bitString'),
    
    out = divVec(x,cut);
    N = de2re(out{1},test.Nmin,test.Nmax);
    lBest = de2re(out{2},test.Lmin,test.Lmax);
    eBest = de2re(out{3},test.Emin,test.Emax);
    
elseif strcmp(type,'doubleVector'),
    N = round(x(1));
    eBest = x(3);
    lBest = x(2);
end

figure(2)
for i=1:size(e,2),
    U=g0*Isp*log((1+l)./(l+e(i)));
    semilogx(l,U)
    hold on
end
U=g0*Isp*log((1+lBest)./(lBest+eBest));
semilogx(lBest,U,'b*')
%legend('toggle')
legend('e = 0.001','e = 0.01','e = 0.05','e = 0.1','e = 0.2','e = 0.5')
 saveas(gcf, strcat(fullfile(PATH, finalName),'.fig'))
end
