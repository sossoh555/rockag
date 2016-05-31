clear all
close all
clc

x = 1:0.1:5;
y = x.^2;
figure('Name','Algo','NumberTitle','off')
plot(x,y)


figure('Name','EITA','NumberTitle','off')
plot(y,x)

fH = get(0,'Children');
f = findobj('name','Algo');
h = get(f);
set(0, 'currentfigure', f);
hold on
plot(2,5,'bo')
hold off
%    plotBest = findobj(get(gca,'Children'),'Tag','gaplotbestf');
%         plotMean = findobj(get(gca,'Children'),'Tag','gaplotmean');
%         newX = [get(plotBest,'Xdata') state.Generation];
%         newY = [get(plotBest,'Ydata') best];
%         set(plotBest,'Xdata',newX, 'Ydata',newY);
%         newY = [get(plotMean,'Ydata') m];
%         set(plotMean,'Xdata',newX, 'Ydata',newY);
%         set(get(gca,'Title'),'String',sprintf('Best: %g Mean: %g',best,m));