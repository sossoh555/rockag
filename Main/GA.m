function GA


clear all
close all
clc
set(0,'DefaultFigureWindowStyle','normal');

global cut Mpay Udes DEBUG g0 Isp type test finalName
global PATH POP fTEST scalingGain fGain
global ALLbool MVECbool DELVbool COSTbool
global W
global Fit P

PATH = 'C:\Users\Bruno\Google Drive\TG\Código\Main\Data\';


test = struct('Nmin',1,'Nmax',6,...
    'Lmin',0.001, 'Lmax',0.1, ...
    'Emin',0.1,'Emax',1);
% fitness

ALLbool = false;
MVECbool = true;
DELVbool = true;
COSTbool = true;
W1 = 1;W2 = 1;W3 = 1;
W = [W1 W2 W3] ;
P1 = 0.5;P2 = 0.5; P3 = 0;
P = [P1 P2 P3];

if sum(P) ~= 1,
   error('P deve ser igual a 1'); 
end

DEBUG = true;


Udes = 5.657;
Isp = 350;
g0 = 9.81/1000;

Fit(3,1) = 0;
Mpay = 10000; %[kg]
Npop = 70; % Tamanho da populacao
Ngen = 100; % Numero de geracoes
Neli = 1; % Numero de elitismo
mutationRate = 0.05; % 5 Percent

nW = 20;
Wbool = true;

if ~Wbool,
    nW = 1;
end

type = 'bitString'; % doubleVector
%type = 'doubleVector';




if strcmp(type,'bitString'),
    var = [5 10 10];
    Nvars = sum(var(:)); %tamanho do indivíduo
    cut = zeros(1,Nvars);
    pos = 0;
    for i = 1:size(var,2)-1,
        pos = pos + var(i);
        cut(pos) = 1;
    end
    
    %[ 5 10 ]
    %[ N Sigma]
    lb = [];
    ub = [];
elseif strcmp(type,'doubleVector')
    Nvars = 3;
    lb = [test.Nmin test.Lmin test.Emin];
    ub = [test.Nmax test.Lmax test.Emax];
else
    error(['Não existe o tipo de GA selecionado: ' type])
end



%=============================
%            FILE
%=============================

DateString = strrep(datestr(datetime('now')),':','');
finalName = strcat(DateString,...
    '-',type,'-',...
    num2str(Npop),'-',...
    num2str(Ngen),'-',...
    num2str(Mpay));


POP = strcat(fullfile(PATH, strcat('POP_',finalName,'.txt')));
fTEST = strcat(fullfile(PATH, strcat('TEST_',finalName,'.txt')));
fGain = strcat(fullfile(PATH, strcat('GAIN_',finalName,'.txt')));
fExcel = strcat(fullfile(PATH, strcat('excelW_',finalName,'.xls')));


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
optsINICIAL = gaoptimset(...
    'PopulationType',type, ...
    'PopulationSize', 70, ...
    'Generations', 2, ...
    'EliteCount', 1, ...
    'CrossoverFcn', @crossoverscattered,...
    'FitnessScalingFcn', @fitscalingrank,...
    'MutationFcn',{@mutationuniform, mutationRate}, ...
    'SelectionFcn', @selectionroulette,...
    'StallGenLimit', 50, ...
    'FitnessScalingFcn',{@fitscalingtop,0.9},...
    'OutputFcn',@outputGA);

scalingGain = true;
Wcheck =[];
for j =1:nW;
    [xbest, fbest, exitflag] = ga(@fitn, Nvars, [], [], [], [], ...
        lb, ub, [], [], optsINICIAL);
    if ~Wbool, 
    fid = fopen(fGain, 'at' );
    fprintf(fid, '\nWFINAL = ');
    fprintf(fid, '%f ',W(:));
    fprintf(fid,'\n');
    end
    
    if Wbool,
        Wcheck = [Wcheck;W];
        disp(W)
        W(:) = 1;
    end 
end
if Wbool,
    for p=1:size(Wcheck,2),
    W(p) = mean(Wcheck(:,p));    
    end
   Wcheck = [Wcheck;W] 
   
  xlswrite(fExcel,Wcheck);
  %a = {1,nW+1,'=sum(a1,a2)'}
  %a(2,:) = {2,nW+1,sprintf('=TEXTO(MÉDIA(B1:B%d);"0,00E+00") & " +-" & TEXTO(DESVPAD.A(B1:B%d);"0,00E+00")',nW)} 
  %a(3,:) = {3,nW+1,sprintf('=TEXTO(MÉDIA(C1:C%d);"0,00E+00") & " +-" & TEXTO(DESVPAD.A(C1:C%d);"0,00E+00")',nW)} 
  %xlswrite(fExcel,a);
end

scalingGain = false;



fid = fopen(POP, 'at' );

fprintf(fid,'\n********* Starting the algorithm *********\n');

ch0 = {'Nvars','Udes', 'Isp', 'Mpay','Npop', 'Ngen', 'Neli','MutatRate'};
ch1 = {'ALLbool', 'DELVbool','MVECbool', 'COSTbool', 'W1', 'W2','W3'};
ch2 = {'Nmin', 'Nmax', 'Lmin', 'Lmax', 'Emin', 'Emax'};

colheadings = [ch0 ch1 ch2];

fms0 = {'d','.2f','d','d','d','d','d','.2f'};
fms1 = {'d','d','d','d','.2E','.2E','.2E'};
fms2 = {'d','d','.2E','.2E','.2E','.2E'};

fms = [fms0 fms1 fms2];

rowheadings = {''};

data0 = [Nvars Udes Isp Mpay Npop Ngen Neli mutationRate];
data1 = [ALLbool DELVbool MVECbool COSTbool W(1) W(2) W(3)];
data2 = [test.Nmin test.Nmax test.Lmin test.Lmax test.Emin test.Emax];
data = [data0 data1 data2];

wid = ones(1,size(colheadings,2));
for i=1:size(colheadings,2),
    wid(i) = length(colheadings{1,i}) ;
end
[truefalse, index] = ismember('Udes', colheadings);
wid(index) = wid(index) + 2;
[truefalse, index] = ismember('W1', colheadings);
wid(index) = wid(index) + 8;
[truefalse, index] = ismember('W2', colheadings);
wid(index) = wid(index) + 8;
[truefalse, index] = ismember('W3', colheadings);
wid(index) = wid(index) + 8;
[truefalse, index] = ismember('Lmin', colheadings);
wid(index) = wid(index) + 2;

displaytable(data,colheadings,wid,fms,rowheadings,fid);

fprintf(fid,'\n\n');




opts = gaoptimset(...
    'PopulationType',type, ...
    'Display','diagnose',...
    'PopulationSize', Npop, ...
    'Generations', Ngen, ...
    'EliteCount', Neli, ...
    'PlotFcns', {@gaplotbestfModified,@gaplotscores,@gaplotdistance,@gaplotscorediversity,@gaplotbestindiv,@gaplotrange},...
    'CrossoverFcn', @crossoverscattered,...
    'FitnessScalingFcn', @fitscalingrank,...
    'MutationFcn',{@mutationuniform, 0.05}, ...
    'SelectionFcn', @selectionroulette,...
    'StallGenLimit', 50, ...
    'FitnessScalingFcn',{@fitscalingtop,0.9},...
    'OutputFcn',@outputGA);

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

%================================
%             GA
%================================
[xbest, fbest, exitflag] = ga(@fitn, Nvars, [], [], [], [], ...
    lb, ub, [], [], opts);
% print(strcat('ALLplots',finalName),'-depsc');
%================================
%           MULTI GA
%================================


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

%analyzeBest(xbest)

%%
% We can also ask |ga| to return the optimal volume of the beam.
%fprintf('\nCost function returned by ga = %g\n', fbest);
%fprintf('%g ', xbest);% fprintf('\n')
% fprintf('%%%%%%%%%%BEST INDIVIDUAL%%%%%%%%%%%\n\n')
% DEBUG = true;
% fitn(xbest);
% DEBUG = false;

%delV = abs(g0*Isp*log((1+xbest(1))/(xbest(2) + xbest(1))));
%fprintf('deltaV: %f \n\n', delV);

%plotN(xbest,exitflag);

fid = fopen(POP, 'at' );
fprintf(fid, 'BEST: \n');
fitn(xbest);

fclose('all');
salvarImagens();
%fclose(fileID);
%analyzeBestVec(xbest)

%% References
%
% [1] Survey of discrete variable optimization for structural design, P.B.
% Thanedar, G.N. Vanderplaats, J. Struct. Eng., 121 (3), 301-306 (1995)
%displayEndOfDemoMessage(mfilename)

end
% function analyzeBestString(xbest)
% global cut Mpay W1 W2 e
%
% g = sprintf('%d ', xbest);
% fprintf('BEST: %s\n', g)
%
% out = divVec(xbest,cut);
% N = de2re(out{1},1,5);
% fprintf('N: %f \n',N)
%
% lambda = de2re(out{2},0.0001,1);
% fprintf('lambda: %f \n',lambda)
%
% obj(1) = abs(real(log(lambda + e*(1 - lambda)))^(-N)/W1);
% fprintf('U: %f \n',obj(1))
%
% Mvec = (Mpay/((lambda)^N));
% fprintf('Vehicle Mass: %f \n',Mvec)
% end

function analyzeBestVec(x)
global W1 e Mpay W2 Udes
g = sprintf('%d ', x);
fprintf('BEST: %s\n', g)

N = x(1);
lambda = x(2);

obj(1) =  abs(Udes - real(log(lambda + e*(1 - lambda))^(N)));

fprintf('U: %f \n',obj(1))

Mvec = (Mpay/((lambda)^N));
fprintf('Vehicle Mass: %f \n',Mvec)
fprintf('Vehicle Mass Pat: %f \n',Mvec/W2)

end

function salvarImagens()
global PATH finalName
h = get(0,'children');

for i=1:length(h)
    name = get(h(i),'name');
    saveas(h(i), strcat(fullfile(PATH, strcat(name,'_',finalName))), 'fig');
end
end

