function FITNESS = fitn(x)
global cut Mpay g0 Isp Udes DEBUG type test

A = 0.6:-0.1:0.1;
B = 15:2:25;
W1 = 1;
W2 = 1E8;
W3 = 4E3;
Cost = 0; 
CostFit = 0;

if strcmp(type,'bitString'),
    
    out = divVec(x,cut);
    N = de2re(out{1},test.Nmin,test.Nmax);
    lambda = de2re(out{2},test.Lmin,test.Lmax);
    e = de2re(out{3},test.Emin,test.Emax);

    if DEBUG,
        fprintf('\n')
        for i =1:size(out,2),
            fprintf('%g ', out{i});
        end
        
    end
    
elseif strcmp(type,'doubleVector'),
    N = round(x(1));
    lambda = x(2);
    e = x(3);
end

%=======================
%          U
%=======================
U = log(((1+lambda)/(e + lambda))^N);
delV = g0*Isp*U;
pPL = lambda/(1+lambda);

num = (pPL^(1/N))*(1-e) + e;
U = log(1/num);
delV = U*Isp*g0*N;

%=======================
%         Mvec
%=======================

mf = Mpay*((lambda+e)/lambda);
num = (pPL^(1/N))*(1-e) + e;
n = (1/num)^N;

Mvec = mf*n;

%Nmax = round(N);
%for i =0:1:Nmax
%mE(i+1) = (1 - pPL^(1/Nmax))*e*Mpay/(pPL^((Nmax - i)/Nmax));
%mp(i+1) = (1 - pPL^(1/Nmax))*(1-e)*Mpay/(pPL^((Nmax - i)/Nmax));
%end
%Mvec = sum(mE)+sum(mp) + Mpay;
%Mvecfit = Mvec/W2;
%Mvec = Mpay*(lambda^(1/N) + 1)/(lambda^(1/N));

%=======================
%        Cost
%=======================
N = round(N);
Cost = 0;
for i=1:N,
Cost = Cost + N*(1 + A(N)*exp(-((0.025*B(N))/(lambda-lambda*e+e))^3.5));
end

%=======================
%        FITNESS
%=======================
delVFit = abs((Udes - delV)/Udes/W1);
MvecFit = Mvec/W2;
CostFit = Cost/W3;
FITNESS = delVFit + MvecFit;% + CostFit;
if FITNESS >1, FITNESS = 1; end

if DEBUG
    %%%%%%%%%%%%%%%% INFORMACOES DO AG %%%%%%%%%%%%%%%%%%%%%%%
    colheadings = {'FIT','DeltaVFit','MvecFit','CostFit','lambda','e','N','DeltaV','Mvec','Cost'};
    rowheadings = {'Res'};
    fms = {'.2E','.2E','.2E','.2E','.4f','.4f','.2f','.4f','.4f','.4f'};
    fileID = 1;
    wid = 9;
    data = [FITNESS delVFit MvecFit CostFit lambda e N delV  Mvec Cost];
    fprintf('\n')
    displaytable(data,colheadings,wid,fms,rowheadings,fileID);
    fprintf('\n')
   
    
    %%%%%%%%%%%%%%%%%% ALGUNS CALCULOS %%%%%%%%%%%%%%%%%%%%%%%%%
    mp = Mvec - mf;
    mE = Mvec - mp - Mpay;
    fprintf('======== CALCULATED ======== \n')
    colheadings = {'N','pPL [%]','Mpay [kg]','mf [kg]','mp [kg]','mE [kg]',...
        'Mvec [kg]','dV [kg/s]','Cost [U$]'};
    rowheadings = {'Res'};
    fms = {'.2f','.4f','.1f','.1f','.1f','.1f','.2f','.4f','.4f'};
    fileID = 1;
    wid = 10;
    data = [N pPL*100 Mpay mf mp mE Mvec delV Cost];
    fprintf('\n')
    displaytable(data,colheadings,wid,fms,rowheadings,fileID);
    fprintf('\n')
   
    
    % fprintf('\n FIT \t Lamb  \t estr \t DelvaV \t Mvec \t Mfit\t\t Cost')
   % fprintf('\n %.3f \t %.3f  \t %.3f \t %.3f \t %.3f \t %.3f \n\n',val, lambda,e,N,delV,Mvec,MvecFit)
end
%
% if DEBUG
%     g = sprintf('%d ', x);
%     fprintf('Answer: %s\n', g)
%     fprintf('N: %f \n',N)
%     fprintf('lambda: %f \n',lambda)
%     fprintf('obj(1): %f \n',obj(1))
%     fprintf('obj(2): %f \n',obj(2))
%     fprintf('val: %f \n',val)
%     fprintf('\n\n')
% end

