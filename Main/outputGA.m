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
global PATH finalName


nameFile =  strcat('POP_',finalName,'.txt');
fid = fopen(strcat(fullfile(PATH, nameFile)), 'at' );

optchanged = false;

pop = state.Population;
scores = state.Score;
best = state.Best;



switch flag
    case 'init'
fn_structdisp(options)
%fprintf(fid,t);
        pop = options.PopulationSize; type = options.PopulationType;
    
    colheadings = {'Type','Pop'};
    rowheadings = {' '};
    fms = {'%s','%d'};
    wid = 10;
    data = [];
    
   % displaytable(data,colheadings,wid,fms,rowheadings,fid);
    
    
        fprintf(fid,'\n********* Starting the algorithm *********\n');
        fprintf(fid,'Population Size: %d\n', options.PopulationSize);
        
        fprintf(fid,'<------> Generation %d <------>\n',state.Generation);
        
        
        for i =1:size(pop,1),
            fprintf(fid,'[%03.f]: ',i);
            fprintf(fid,'%g ',pop(i,:));
            fprintf(fid,'%f',scores(i,1));
            fprintf(fid,'\n');
        end
    case {'iter','interrupt'}
        
        fprintf(fid,'Iterating ...');
        fprintf(fid,'<------> Generation %d <------>\n',state.Generation);
        for i =1:size(pop,1),
            fprintf(fid,'[%03.f]: ',i);
            fprintf(fid,'%g ',pop(i,:));
            fprintf(fid,'%f',scores(i,1));
            if scores(i,1) == best(end),  fprintf(fid,' [BEST] '); end
            fprintf(fid,'\n');
        end
    case 'done'
        fprintf(fid,'Performing final task');
end
fclose(fid);
end
