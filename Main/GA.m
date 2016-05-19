clear all 
close all
clc

global cut Mpay 

N = 15; % Tamanho do Individuo
Npop = 50; % Tamanho da populacao


Mpay = 100; %[kg]
cut = zeros(1,N);
cut(5) = 1; % corte do individuo
%[ 5 10 ]
%[ N Sigma]

%lb = [1 30 2.4 45 2.4 45 1 30 1 30];
%ub = [5 65 3.1 60 3.1 60 5 65 5 65];

%%
%
% _Set the Options_
%
% To obtain a more accurate solution, we increase the PopulationSize, and
% Generations options from their default values, and decrease the
% EliteCount and TolFun options. These settings cause |ga| to use a larger
% population (increased PopulationSize), to increase the search of the
% design space (reduced EliteCount), and to keep going until its best
% member changes by very little (small TolFun).  We also specify a plot
% function to monitor the penalty function value as |ga| progresses.
%
% Note that there are a restricted set of |ga| options available when
% solving mixed integer problems - see Global Optimization Toolbox User's
% Guide for more details.
opts = gaoptimset(...
    'PopulationType','bitstring', ...   
    'PopulationSize', Npop, ...
    'Generations', 50, ...
    'EliteCount', 1, ...%'TolFun', 1e-8, ...
    'PlotFcns', @gaplotbestfModified);

%%
%
% _Call |ga| to Solve the Problem_
%
% We can now call |ga| to solve the problem. In the problem statement $x_1$
% and $x_2$ are integer variables. We specify this by passing the index
% vector |[1 2]| to |ga| after the nonlinear constraint input and before
% the options input. We also seed and set the random number generator here
% for reproducibility.
rng(0, 'twister');
[xbest, fbest, exitflag] = gamultiobj(@fitn, N, [], [], [], [], ...
    [], [], [], [], opts);
    
%% 
%
% _Analyze the Results_
%
% If a problem has integer constraints, |ga| reformulates it internally. In
% particular, the fitness function in the problem is replaced by a penalty
% function which handles the constraints. For feasible population members,
% the penalty function is the same as the fitness function.
%
% The solution returned from |ga| is displayed below. Note that the section
% nearest the support is constrained to have a width ($x_1$) and height
% ($x_2$) which is an integer value and this constraint has been honored by
% GA.
display(xbest);

%%
% We can also ask |ga| to return the optimal volume of the beam. 
fprintf('\nCost function returned by ga = %g\n', fbest);


%% References
%
% [1] Survey of discrete variable optimization for structural design, P.B.
% Thanedar, G.N. Vanderplaats, J. Struct. Eng., 121 (3), 301-306 (1995)
displayEndOfDemoMessage(mfilename)