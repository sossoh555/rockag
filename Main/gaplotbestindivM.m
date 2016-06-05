function state = gaplotbestindivM(options,state,flag)
%GAPLOTBESTINDIV Plots the best individual.
%   STATE = GAPLOTBESTINDIV(OPTIONS,STATE,FLAG) plots the best 
%   individual's genome as a histogram, with the number of bins
%   in the histogram equal to the length of the genome.
%
%   Example:
%    Create an options structure that uses GAPLOTBESTINDIV
%    as the plot function
%     options = gaoptimset('PlotFcns',@gaplotbestindiv);

%   Copyright 2003-2007 The MathWorks, Inc.

if  size(state.Score,2) > 1
    title('Gráfico do Melhor Indivíduo: não disponível','interp','none');
    return;
end

switch flag
    case 'init'
        GenomeLength = size(state.Population,2);
        [~,i] = min(state.Score);
        h = bar(double(state.Population(i,:)));
        set(h,'edgecolor','none','Tag','gaplotbestindiv')
        set(gca,'xlim',[0,1 + GenomeLength])
        title('Melhor Indivíduo Atual','interp','none')
        xlabel(sprintf('Número de Variáveis (%i)',GenomeLength),'interp','none');
        ylabel('Valor do Gene','interp','none');

    case 'iter'
        [~,i] = min(state.Score);
        h = findobj(get(gca,'Children'),'Tag','gaplotbestindiv');
        set(h,'Ydata',double(state.Population(i,:)));
end


