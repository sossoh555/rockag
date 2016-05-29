clear all
close all
clc

opt = 'ALL'
%LAST
opengl('save','software')
N = 1:1:15;
pPL = 0.05;
e = 0.15;
Isp = 350;
g = 9.81/1000;

num = (pPL.^(1./N))*(1-e) + e;
vb = log(1./num)*Isp*g.*N;

hPlot = plot(N,vb,'--bs',...
    'LineWidth',1.2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1,1,1])
hXLabel = xlabel('N')
hYLabel = ylabel('\DeltaV (km/s)');
hTitle  = title ('Escrever Título');

vbinf = Isp*g*(1-e)*log(1/pPL);
text(N(end)-2, vbinf,strcat('v_{b\infty} =', {'  '}, num2str(vbinf), {' km/s'}),...
    'VerticalAlignment','bottom')

vbinf(1:1:size(N,2)+1) = vbinf;
hLimit = line([0 N],vbinf);

set(hLimit,'LineWidth',2,'Color','r');


if opt == 'ALL',
    
    strValues = strtrim(cellstr(num2str([vb(:)],'%.2f')));
    hText = text(N,vb,strValues,'HorizontalAlignment','right','VerticalAlignment','bottom');
    
elseif opt == 'LAST'
    strValues = strtrim(cellstr(num2str([vb(end)],'%.7f')));
    hText = text(N(end),vb(end),strValues,'HorizontalAlignment','right','VerticalAlignment','bottom');
end



set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set(hText,'FontName','AvantGarde','FontSize',11);
set([hXLabel, hYLabel],'FontSize',10);
set(hTitle,'FontSize', 12,'FontWeight','bold');
set(gca,'XTick', 0:1:N(end)+ 1);

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:1:10, ...
    'LineWidth'   , 1         );
grid on
