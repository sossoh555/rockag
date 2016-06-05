function state = gaplotdistanceM(options,state,flag)
%GAPLOTDISTANCE Averages several samples of distances between individuals.
%   STATE = GAPLOTDISTANCE(OPTIONS,STATE,FLAG) plots an averaged distance
%   between individuals.
%
%   Example:
%    Create an options structure that uses GAPLOTDISTANCE
%    as the plot function
%     options = gaoptimset('PlotFcns',@gaplotdistance);

%   Copyright 2003-2007 The MathWorks, Inc.

samples = 20;
choices = ceil(sum(options.PopulationSize) * rand(samples,2));
switch flag
    case 'init'
        population = state.Population;
        distance = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            distance = distance + sqrt( sum ( d.* d));
        end
        plotDist = plot(state.Generation,distance/samples,'.');
        set(gca,'xlimmode','manual','zlimmode','manual', ...
            'alimmode','manual')
        set(gca,'xlim',[1,options.Generations]);
        set(plotDist,'Tag','gaplotdistance');
        xlabel('Geração','interp','none');
        ylabel('Distância Média');
        title('Distância média entre os Indivíduos','interp','none')

    case 'iter'
        population = state.Population;
        distance = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            distance = distance + sqrt( sum ( d.* d));
        end
        plotDist = findobj(get(gca,'Children'),'Tag','gaplotdistance');
        newX = [get(plotDist,'Xdata') state.Generation];
        newY = [get(plotDist,'Ydata') distance/samples];
        set(plotDist,'Xdata',newX,'Ydata',newY);
end
