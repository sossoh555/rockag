clear all
close all
clc

x = 1:0.1:5;
y = x.^2;
figure('Name','Algo','NumberTitle','off','GraphicsSmoothing', 'on')
dataPlotN = plot(x,y);


for i=1:10,

%figure('Name','EITA','NumberTitle','off')
%h('Tag', 'Eita') = semilogx(y,x);

%fH = get(0,'Children');

f = findobj('name','Algo');
set(0, 'currentfigure', f);
hold on
h(i) = plot(3,i,'bx','Color',[0,0,0])



 %newX = [get(dataPlotN,'Xdata') 1];
 %newY = [get(dataPlotN,'Ydata') 4];
 %semilogx(newX,newY)
end
g = findobj('type','line')
 %set(dataPlotN,'Xdata',newX, 'Ydata',newY);

%    best = min(state.Score);
%         m    = meanf(state.Score);
%         %set(gca,'xlim',[0,state.Generation + 5]);
%         plotBest = findobj(get(gca,'Children'),'Tag','gaplotbestf');
%         plotMean = findobj(get(gca,'Children'),'Tag','gaplotmean');
%         newX = [get(plotBest,'Xdata') state.Generation];
%         newY = [get(plotBest,'Ydata') best];
%         set(plotBest,'Xdata',newX, 'Ydata',newY);
%         newY = [get(plotMean,'Ydata') m];
%         set(plotMean,'Xdata',newX, 'Ydata',newY);
%         set(get(gca,'Title'),'String',sprintf('Best: %g Mean: %g',best,m));;