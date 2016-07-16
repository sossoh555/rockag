function [state, options,optchanged] = outputGA(options,state,flag)

%GAOUTPUTFCNTEMPLATE Template to write custom OutputFcn for GA.
%   [STATE, OPTIONS, OPTCHANGED] = GAOUTPUTFCNTEMPLATE(OPTIONS,STATE,FLAG)
%   where OPTIONS is an options structure used by GA.
%
%   STATE: A structure containing the following information about the state
%   of the optimization:
%             Population: Population in the current generation
%                  Score: Scores of the current population
%             Generation: Current generation number
%              StartTime: Time when GA started
%               StopFlag: String containing the reason for stopping
%              Selection: Indices of individuals selected for elite,
%                         crossover and mutation
%            Expectation: Expectation for selection of individuals
%                   Best: Vector containing the best score in each generation
%        LastImprovement: Generation at which the last improvement in
%                         fitness value occurred
%    LastImprovementTime: Time at which last improvement occurred
%
%   FLAG: Current state in which OutputFcn is called. Possible values are:
%         init: initialization state
%         iter: iteration state
%    interrupt: intermediate state
%         done: final state
%
%   STATE: Structure containing information about the state of the
%          optimization.
%
%   OPTCHANGED: Boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, GA, GAOPTIMSET

%   Copyright 2004-2006 The MathWorks, Inc.
global POP fTEST
global Fit P W  scalingGain fGain

if scalingGain,
    fid = fopen(fGain, 'at' );
else
    fid = fopen(POP, 'at' );
end

optchanged = false;

pop = state.Population;
scores = state.Score;
best = state.Best;

switch flag
    case 'init'
        if ~scalingGain,
            [~,i] = min(state.Score);
            plotN(state.Population(i,:),flag);
        end
        Fit(:) = 0;
        
        pop = options.PopulationSize;
        typeGA = options.PopulationType;
        gen = options.Generations;
        mutat = func2str(options.MutationFcn);
        creat = func2str(options.CreationFcn);
        fitnessScaling = func2str(options.FitnessScalingFcn);
        selec = func2str(options.SelectionFcn);
        cross = func2str(options.CrossoverFcn);
        tolFun = options.TolFun;
        tolCon = options.TolCon;
        crossFraction = options.CrossoverFraction;
        elite = options.EliteCount;
        migDir = options.MigrationDirection;
        migInt = options.MigrationInterval;
        migFrac =options.MigrationFraction;
        fncs = options.PlotFcns;
        plotFcns = '';
        for i =1:size(fncs,2),
            plotFcns =  [plotFcns  sprintf(' %s ',func2str(fncs{i}))];
        end
        names = {'Type' 'Population Size' 'Mutation' 'Creation' ...
            'Fitness Scaling' 'Selection' 'Crossover' 'PlotFunctions'...
            'TolFun' 'TolCon' 'Crossover Fraction' 'Elite Count'...
            'Migration Direction' 'Migration Interval'...
            'Migration Fraction'};
        maxChars = 0;
        for i=1:size(names,2)
            ref = length(names{i});
            if maxChars < ref, maxChars = ref;   end
        end
        
        fprintf(fid,'Type: %s\n', typeGA);
        fprintf(fid,'Population Size: %d\n',pop);
        fprintf(fid,'Generations: %d\n',gen);
        fprintf(fid,'Mutation: %s\n',mutat);
        fprintf(fid,'Creation: %s\n',creat);
        fprintf(fid,'Fitness Scaling: %s\n',fitnessScaling);
        fprintf(fid,'Selection: %s\n',selec);
        fprintf(fid,'Crossover: %s\n',cross);
        fprintf(fid,'PlotFunctions: %s\n',plotFcns);
        fprintf(fid,'Crossover Fraction: %.2f\n',crossFraction);
        fprintf(fid,'TolFun: %.2E\n',tolFun);
        fprintf(fid,'TolCon: %.2E\n',tolCon);
        fprintf(fid,'Elite Count: %d\n',elite);
        fprintf(fid,'Migration Direction: %s\n',migDir);
        fprintf(fid,'Migration Interval: %d\n',migInt);
        fprintf(fid,'Migration Fraction: %.3f\n\n\n',migFrac);
        
        clearvars -global Fit;
        Fit(3,1) = 0;
        
    case {'iter','interrupt'}
        
        widp = 3;
        fprintf(fid,'Iterating ...');
        fprintf(fid,'<------> Generation %03d <------>\n',state.Generation);
        
        div = zeros(1,size(pop,2));
        
        for i =1:size(pop,1),
            fprintf(fid,'[%03.f]: ',i);
            for j =1:size(pop,2),
                fprintf(fid,'%*g ',widp,pop(i,j));
                if pop(i,j) == 1,
                    div(j) = div(j) + 1;
                end
            end
            
            fprintf(fid,'| %f',scores(i,1));
            if scores(i,1) == best(end),
                fprintf(fid,' [BEST] |');
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'[DIV]: ');
        for k =1:size(div,2),
            aux = div(1,k);
            
            fprintf(fid,'%*g ',widp,aux);
        end
        
        if ~scalingGain && state.LastImprovement == state.Generation,
            [~,i] = min(state.Score);
            plotN(state.Population(i,:),flag);
            
            fprintf(fid,'\n');
            fprintf(fid, 'BEST Geracao %d:\n',state.Generation);
            fitn(state.Population(i,:));
            fid = fopen(POP, 'at' );
            fprintf(fid,'\n');
            fprintf(fid,'\n');
            
        end
        
        if scalingGain,
            fid = fopen(fGain, 'at' );
        else
            fid = fopen(fTEST, 'at' );
        end
        FitAvg = zeros(1,3);
        
        if state.Generation == 1 && scalingGain,
            
            FitAvg(1) = sum(Fit(1,1:end))/options.PopulationSize;
            FitAvg(2) = sum(Fit(2,1:end))/options.PopulationSize;
            FitAvg(3) = sum(Fit(3,1:end))/options.PopulationSize;
            for i=1:size(FitAvg,2),
                fprintf(fid, 'Fit%d: %f\n',i,FitAvg(i));
            end
            
            C = sum(FitAvg);
            fprintf(fid, 'C: %f\n',C);
            count = 0;
            aux = 0;
            for i=1:size(FitAvg,2),
                fprintf(fid, 'W%d_inicial: %f\n',i,W(i));
                for j=1:count,
                    aux = aux + FitAvg(j)*W(j);
                end
                fprintf(fid, 'aux: %f\n',aux);
                W(i) = (sum(P(1:(count+1)))*C - aux)/FitAvg(i);
                fprintf(fid, 'W%d_final: %f\n',i,W(i));
                count = count + 1;
                aux = 0;
            end
            
            fprintf(fid, 'Geracao: %d\n',state.Generation);
            for i=1:size(W,2),
                fprintf(fid,'W%d: %f\n',i,W(i));
            end
            fprintf(fid,'\n');
        else
            fprintf(fid, 'Geracao: %d\n',state.Generation);
            
            for i=1:size(FitAvg,2),
                fprintf(fid,'Fit #%d\n',i);
                
                fprintf(fid, 'Min: %.3E\n',min(Fit(i,2:end)));
                fprintf(fid, 'Max: %.3E\n',max(Fit(i,2:end)));
            end
        end
        
        clearvars -global Fit;
        Fit(3,1) = 0;
    case 'done'
        if ~scalingGain,
            [~,i] = min(state.Score);
            plotN(state.Population(i,:),flag);
        end
        fprintf(fid,'\nPerforming final task\n');
end
fclose('all');
end
