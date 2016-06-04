clear all
close all
clc

x = 1:0.1:5;
y = x.^2;
z = x+4;
figure('Name','Algo','NumberTitle','off')
plot(x,y); hold on
plot(x,z,'bx');

mx = findobj(get(gca,'Children'),'Marker','x')

for i=1:10,

%figure('Name','EITA','NumberTitle','off')
%h('Tag', 'Eita') = semilogx(y,x);

%fH = get(0,'Children');

f = findobj('name','Algo');
set(0, 'currentfigure', f);



hold on
%h(i) = plot(3,i,'bx','Color',[0,0,0]);


 %newX = [get(dataPlotN,'Xdata') 1];
 %newY = [get(dataPlotN,'Ydata') 4];
 %semilogx(newX,newY)
end

axesO = get(gca,'children')
newX = [get(mx,'Xdata') 7];
newY = [get(mx, 'Ydata') 30];
set(mx, 'Xdata',newX, 'Ydata', newY)

if axesO(1).Marker == '-',
    disp('EH PORRA')
else
    disp('NAO EH PORRA')
end


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