function FITNESS = fitnMulti(x)
global cut Mpay g0 Isp Udes DEBUG type test POP
global ALLbool MVECbool DELVbool COSTbool scalingGain fGain
global Fit

A = 0.6:-0.1:0.1;
B = 15:2:25;
Cost = 0;delV = 0;Mvec = 0; mf = 0; mp = 0; mE = 0;
CostFit = 0;delVFit = 0;MvecFit = 0;

if strcmp(type,'bitString'),

    out = divVec(x,cut);
    N = de2re(out{1},test.Nmin,test.Nmax);
    lambda = de2re(out{2},test.Lmin,test.Lmax);
    %lambda = lambda^N;
    e = de2re(out{3},test.Emin,test.Emax);
    
%     if DEBUG,
%         fprintf('\n')
%         for i =1:size(out,2),
%             fprintf('%g ', out{i});
%         end
%         
%     end
    
elseif strcmp(type,'doubleVector'),
    N = round(x(1));
    lambda = x(2);
    e = x(3);
end

%=======================
%Important Considerations
%=======================
pPL = (lambda/(1+lambda))^N;


%=======================
%          delV
%=======================

if DELVbool || ALLbool,
    % U = log(((1+lambda)/(e + lambda))^N);
    %delV = g0*Isp*U;

    
    num = (pPL^(1/N))*(1-e) + e;
    U = log(1/num);
    delV = U*Isp*g0*N;
    %num = e*((1/pPL)^(1/N)) + (1-e);
    %delV = Isp*g0*N*((1/pPL)^(1/N))/num;
    
    delVFit = abs((Udes - delV)/Udes);
end




%=======================
%         Mvec
%=======================

if MVECbool || ALLbool,
%     mf = Mpay*((lambda+e)/lambda);
%     num = (pPL^(1/N))*(1-e) + e;
%     n = (1/num)^N;
%     Mvec = mf*n;
%     MvecFit = Mvec*W(2);


Nmax = round(N);
for i =0:1:Nmax
mE(i+1) = (1 - pPL^(1/Nmax))*e*Mpay/(pPL^((Nmax - i)/Nmax));
mp(i+1) = (1 - pPL^(1/Nmax))*(1-e)*Mpay/(pPL^((Nmax - i)/Nmax));
end
Mvec = sum(mE)+sum(mp) + Mpay;
MvecFit = Mvec;
%Mvec = Mpay*(lambda^(1/N) + 1)/(lambda^(1/N));

end



%=======================
%        Cost
%=======================
if COSTbool || ALLbool,
    Nr = round(N);
    Cost = 0;
    
    for i=1:Nr,
        A = -0.1*Nr + 0.7;
        B = 2*Nr + 13;
        L = lambda^Nr;
        Cost = Cost + Nr*(1 + A*exp(-((0.025*B)/(L-L*e+e))^3.5));
    end
    CostFit = Cost;
end


%=======================
%        FITNESS
%=======================
% delVFit = delVFit*100;
% MvecFit = MvecFit*100;
% CostFit = CostFit*100;
dir = size(Fit,2);
%fprintf('Dir: %d\n', dir)

Fit(1,dir+1) = delVFit;% [r c] = size(Fit); fprintf('Size1: [%d %d]\n', r,c )
Fit(2,dir+1) = MvecFit;% [r c] = size(Fit); fprintf('Size1: [%d %d]\n', r,c )
Fit(3,dir+1) = CostFit;% [r c] = size(Fit); fprintf('Size1: [%d %d]\n', r,c )


FITNESS(1) = delVFit;
FITNESS(2) = MvecFit;
FITNESS(3) = CostFit;

%FITNESS = delVFit + MvecFit + CostFit;
%if FITNESS >1, FITNESS = 1; end

if DEBUG,
    %nameFile =  strcat(POP,finalName,'.txt');
    %fid = fopen(strcat(fullfile(PATH, nameFile)), 'at' );
    if scalingGain,
        fid = fopen(fGain, 'at' );
    else
        fid = fopen(POP, 'at' );
    end
    
    fprintf(fid,'%g ',x(:));
    fprintf(fid,'| %f',FITNESS);
    fprintf(fid,'\n');
    %%%%%%%%%%%%%%%% INFORMACOES DO AG %%%%%%%%%%%%%%%%%%%%%%%
    mp = sum(mp);
    mE = sum(mE);
    mf = Mvec - mp;
    %mp = Mvec - mf;
    %mE = Mvec - mp - Mpay;
    
    ch1 = {'FIT','dVFit','MvecFit','CFit','dVF [%]','MvF [%]','CF [%]','N','lambda','e','delV [kg/s]','Mvec [kg]','Cost [U$]'};
    ch2 = {'pPL [%]','Mpay [kg]','mf [kg]','mp [kg]','mE [kg]'};
    
    colheadings = [ch1 ch2];
    
    fms1 = {'.2E','.2E','.2E','.2E','.3f','.3f','.3f','.2f','.4f','.4f','.4f','.2E','.4f'};
    fms2 = {'.4f','.2E','.2E','.2E','.2E'};
    
    fms = [fms1 fms2];
    
    rowheadings = {''};
    
    wid = ones(1,size(colheadings,2));
    for i=1:size(colheadings,2),
        wid(i) = length(colheadings{1,i}) ;
    end
    
    [truefalse, index] = ismember('MvecFit', colheadings);
    wid(index) = wid(index) + 3;
    [truefalse, aux] = ismember('FIT', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('dVFit', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('CFit', colheadings);wid(aux) = wid(index);
    
    [truefalse, index] = ismember('N', colheadings);
    wid(index) = wid(index)+3;
    
    [truefalse, index] = ismember('lambda', colheadings);
    [truefalse, aux] = ismember('e', colheadings);wid(aux) = wid(index);
    
    
    
    wid(end) = wid(end) + 1;
    wid(end-1) = wid(end-1) + 1;
    wid(end-2) = wid(end-2) + 1;
    wid(end-3) = wid(end-3) + 1;
    somaFit = sum(FITNESS);
    
    dVFp = 100*delVFit/somaFit;
    MvFp = 100*MvecFit/somaFit;
    CFp = 100*CostFit/somaFit;
    
    
    data1 = [somaFit delVFit MvecFit CostFit dVFp MvFp CFp N lambda e delV Mvec Cost];
    data2 = [pPL*100 Mpay mf mp mE];
    data = [data1 data2];
    
    
    displaytable(data,colheadings,wid,fms,rowheadings,fid);
    fprintf(fid,'\n');
    fclose('all');
    
    
    
    
    
    
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

