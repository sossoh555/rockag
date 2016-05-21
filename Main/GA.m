function GA

clear all 
close all
clc

global  cut Mpay Udes DEBUG

%lobal cut Mpay Udes

var = [5 8];
Nvars = sum(var(:)); % Tamanho do Individuo
Npop = 10; % Tamanho da populacao
Ngen = 10; % Numero de geracoes
Neli = 1; % Numero de elitismo

DEBUG = false;
Udes = 5; 
Mpay = 5000; %[kg]
cut = zeros(1,Nvars);
cut(var(1:end-1)) = 1; % corte do individuo
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
    'Display','diagnose',...
    'PopulationSize', Npop, ...
    'Generations', Ngen, ...
    'EliteCount', Neli, ...
    'PlotFcns', {@gaplotbestfModified,@gaplotscores,@gaplotdistance},...
    'MutationFcn',{@mutationuniform,0.05}, ...
    'SelectionFcn', @selectionroulette);

%%
%
% _Call |ga| to Solve the Problem_
%
% We can now call |ga| to solve the problem. In the problem statement $x_1$
% and $x_2$ are integer variables. We specify this by passing the index
% vector |[1 2]| to |ga| after the nonlinear constraint input and before
% the options input. We also seed and set the random number generator here
% for reproducibility.
%rng(0, 'twister');
[xbest, fbest, exitflag] = ga(@fitn, Nvars, [], [], [], [], ...
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
%display(xbest);

analyzeBest(xbest)

%%
% We can also ask |ga| to return the optimal volume of the beam. 
fprintf('\nCost function returned by ga = %g\n', fbest);


%% References
%
% [1] Survey of discrete variable optimization for structural design, P.B.
% Thanedar, G.N. Vanderplaats, J. Struct. Eng., 121 (3), 301-306 (1995)
%displayEndOfDemoMessage(mfilename)

end 
function analyzeBest(xbest)
global cut Mpay W1 W2 e

g = sprintf('%d ', xbest);
fprintf('BEST: %s\n', g)

out = divVec(xbest,cut);
N = de2re(out{1},1,5);
fprintf('N: %f \n',N)

lambda = de2re(out{2},0.0001,1);
fprintf('lambda: %f \n',lambda)

obj(1) = abs(real(log(lambda + e*(1 - lambda)))^(-N)/W1);
fprintf('U: %f \n',obj(1))

Mvec = (Mpay/((lambda)^N));
fprintf('Vehicle Mass: %f \n',Mvec)


end 
