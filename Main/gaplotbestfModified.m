function state = gaplotbestfModified(options,state,flag)
%GAPLOTBESTF Plots the best score and the mean score.
%   STATE = GAPLOTBESTF(OPTIONS,STATE,FLAG) plots the best score as well
%   as the mean of the scores.
%
%   Example:
%    Create an options structure that will use GAPLOTBESTF
%    as the plot function
%     options = gaoptimset('PlotFcns',@gaplotbestf);

%   Copyright 2003-2007 The MathWorks, Inc.

if size(state.Score,2) > 1
    title('Best Fitness Plot: not available','interp','none');
    return;
end

switch flag
    case 'init'
        hold on;
        %set(gca,'xlim',[0,state.Generation + 5]);
        
        xlabel('Geração','interp','none');
        ylabel('Valor do Fitness','interp','none');
        plotBest = plot(state.Generation,min(state.Score),'.k');
        set(plotBest,'Tag','gaplotbestf');
        %%plotMean = plot(state.Generation,meanf(state.Score),'.b');
        %%set(plotMean,'Tag','gaplotmean');
        plotChange = plot(state.Generation, min(state.Score),'*r','MarkerSize',5);
        set(plotChange, 'Tag','gaplotchange');
        title(['Melhor: ',' Média: '],'interp','none')
    case 'iter'
        best = min(state.Score);
        m    = meanf(state.Score);
        %set(gca,'xlim',[0,state.Generation + 5]);
        plotBest = findobj(get(gca,'Children'),'Tag','gaplotbestf');
        %%plotMean = findobj(get(gca,'Children'),'Tag','gaplotmean');
        plotChange = findobj(get(gca,'Children'), 'Tag', 'gaplotchange');
        newX = [get(plotBest,'Xdata') state.Generation];
        newY = [get(plotBest,'Ydata') best];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
        if best < newY(end-1),
            newXchange = [get(plotChange,'Xdata') state.Generation];
            newYchange = [get(plotChange,'Ydata') best];
            set(plotChange, 'Xdata',newXchange, 'Ydata',newYchange);
        else
            
            %set(plotBest,'Xdata',newX, 'Ydata',newY, 'Color', 'k');
        end
        
        %%newY = [get(plotMean,'Ydata') m];
       %%set(plotMean,'Xdata',newX, 'Ydata',newY);
        set(get(gca,'Title'),'String',sprintf('Melhor: %g Média: %g',best,m));
        
        
    case 'done'
        %LegnD = legend('Melhor Fitness', 'Melhor indivíduo Modificado');
        LegnD = legend('Melhor Fitness', 'Melhor indivíduo Modificado');
        set(LegnD,'FontSize',8);
        hold off;
end

%------------------------------------------------
function m = meanf(x)
nans = isnan(x);
x(nans) = 0;
n = sum(~nans);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(x) ./ n;

