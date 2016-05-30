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
global POP


%nameFile =  strcat(POP,finalName,'.txt');
%fid = fopen(strcat(fullfile(PATH, nameFile)), 'at' );
fid = fopen(POP, 'at' );

optchanged = false;

pop = state.Population;
scores = state.Score;
best = state.Best;



switch flag
    case 'init'
        %fprintf(fid,t);
        
        
        %fprintf(fid, 'População Inicial\n')
        %fprintf(fid,'<------> Generation %03d <------>\n',state.Generation);
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
        
        
        names = {'Type' 'Population Size' 'Mutation' 'Creation' 'Fitness Scaling' 'Selection' 'Crossover' 'PlotFunctions' 'TolFun' 'TolCon' 'Crossover Fraction' 'Elite Count' 'Migration Direction' 'Migration Interval' 'Migration Fraction'};
        maxChars = 0;
        for i=1:size(names,2)
            ref = length(names{i});
            if maxChars < ref, maxChars = ref;   end
        end
        
        %data = {typeGA pop mutat creat fitnessScaling selec cross tolFun tolCon crossFraction elite migDir migInt migFrac};
        
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
        
        
        
        %fprintf(fid,'<------> Generation %03d <------>\n',state.Generation);
        
        
        %         for i =1:size(pop,1),
        %             fprintf(fid,'[%03.f]: ',i);
        %             fprintf(fid,'%g ',pop(i,:));
        %             fprintf(fid,'%f',scores(i,1));
        %             fprintf(fid,'\n');
        %         end
        
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
            if scores(i,1) == best(end),  fprintf(fid,' [BEST] |'); end
            fprintf(fid,'\n');
        end
            fprintf(fid,'[XXX]: ');
            for k =1:size(div,2),
              aux = div(1,k);
              
              fprintf(fid,'%*g ',widp,aux); 
            end
            fprintf(fid,'\n');
        fprintf(fid,'\n');
    case 'done'
        fprintf(fid,'Performing final task');
end
fclose('all');
end
