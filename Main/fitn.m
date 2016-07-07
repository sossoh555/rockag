function FITNESS = fitn(x)
global cut Mpay g0 Isp Udes DEBUG type test POP
global ALLbool MVECbool DELVbool COSTbool scalingGain fGain
global W
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
 
elseif strcmp(type,'doubleVector'),
    N = round(x(1));
    lambda = x(2);
    e = x(3);
end

%=======================
%Important Considerations
%=======================
%pPL = (lambda/(1+lambda));

%=======================
%          delV
%=======================

if DELVbool || ALLbool,
    [pPL, delV] = delVcalc(lambda,e,N);
    delVFit = abs(W(1)*(Udes - delV)/Udes);
end

%=======================
%         Mvec
%=======================
%if MVECbool || ALLbool,

    Nmax = round(N);
    for i =0:1:Nmax - 1
        pPL = (lambda/(1+lambda))^Nmax;
        mE(i+1) = (1 - pPL^(1/Nmax))*e*Mpay/(pPL^((Nmax - i)/Nmax));
        mp(i+1) = (1 - pPL^(1/Nmax))*(1-e)*Mpay/(pPL^((Nmax - i)/Nmax));
    end
    Mvec = sum(mE)+sum(mp) + Mpay;

    MvecFit = Mvec*W(2);
    
    
%end
%=======================
%        Cost
%=======================
%if COSTbool || ALLbool,
    Nr = floor(N);
    rem = mod(N,Nr);
    LOX = 0.08;
    RP = 0.20;
    costE = [0.6 0.2 0.1 0.05 0.025 0.025];
    RPLOX = 0.430585683;
    
    Nr = round(N);
    Cost = 0;
    %
    %     for i=1:Nr,
    %         A = -0.1*Nr + 0.7;
    %         B = 2*Nr + 13;
    %         L = lambda^Nr;
    %         Cost = Cost + Nr*(1 + A*exp(-((0.025*B)/(L-L*e+e))^3.5));
    %     end
    %     CostFit = Cost;
    % end
    for i=1:Nr,
        A = -0.1*Nr + 0.7;
        B = 2*Nr + 13;
%         L = lambda/(lambda + 1);
%         E = e/(e + 1);
%         L = lambda;
%         E = e;
        mi = mE(Nr);
        %Cost = Cost + 1 + A*exp(-mi/400000);
         
        Cost = Cost + mp(Nr)*RPLOX*RP + mp(Nr)*(1-RPLOX)*LOX;
        %Cost = Cost + Nr*(1 + A*exp(-((0.025*B)/(L-L*E+E))^3.5));
    end
    A = -0.1*N + 0.7;
%     B = 2*N + 13;
%     L = lambda/(lambda + 1);
%     E = e/(e + 1);
%     L = lambda;
%     E = e;
    mi = mE(Nr)*rem;
    %Cost = Cost + (1 + A*exp(-mi/400000))*rem;
     Cost = Cost + (mp(Nr)*RPLOX*RP + mp(Nr)*(1-RPLOX)*LOX)*rem;
    Cost = Cost + 299 + 14.2*log(sum(mE));
    for j = 1:Nr,
    
    end
     CostFit = Cost*W(3);
%end


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

% if ~MVECbool && COSTbool,
%     MvecFit = 0;
% end
if ~MVECbool, MvecFit = 0; end
if ~COSTbool, CostFit = 0; end
if ~DELVbool, delVFit = 0; end
% MvecFit = 0;
fitness = delVFit + MvecFit + CostFit;
FITNESS = 100 + log10(fitness);
%if FITNESS > 1, FITNESS = 1; end


Penalty = 0;
fPenal = Penalty*pPL/0.05;
FITNESS = FITNESS + fPenal;
if DEBUG,
    %nameFile =  strcat(POP,finalName,'.txt');
    %fid = fopen(strcat(fullfile(PATH, nameFile)), 'at' );
    if scalingGain,
        fid = fopen(fGain, 'at' );
    else
        fid = fopen(POP, 'at' );
    end
    
    fprintf(fid,'%g ',x(:));
    fprintf(fid,'| %f |%.2E',FITNESS,fitness);
    fprintf(fid,'\n');
    %%%%%%%%%%%%%%%% INFORMACOES DO AG %%%%%%%%%%%%%%%%%%%%%%%
    mp = sum(mp);
    mE = sum(mE);
    mf = Mvec - mp;
    %mp = Mvec - mf;
    %mE = Mvec - mp - Mpay;
    
    ch1 = {'FIT','Penal','dVFit','MvecFit','CFit','dVF [%]','MvF [%]',...
        'CF [%]','N','lambda','e','delV [kg/s]','Mvec [kg]','Cost [U$]'};
    ch2 = {'pPL [%]','Mpay [kg]','mf [kg]','mp [kg]','mE [kg]'};

    colheadings = [ch1 ch2];
    
    fms1 = {'.2E','.2E','.2E','.2E','.2E','.3f','.3f','.3f','.2f','.4f',...
        '.4f','.7f','.3f','.4f'};
    fms2 = {'.4f','.2E','.5E','.5E','.5E'};
    fms = [fms1 fms2];
    
    rowheadings = {''};
    
    wid = ones(1,size(colheadings,2));
    for i=1:size(colheadings,2),
        wid(i) = length(colheadings{1,i}) ;
    end
    
    [truefalse, index] = ismember('mE [kg]', colheadings);
    wid(index) = wid(index) + 5;
    [truefalse, aux] = ismember('mp [kg]', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('mf [kg]', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('Mvec [kg]', colheadings);wid(aux) = wid(index);
    
    [truefalse, index] = ismember('MvecFit', colheadings);
    wid(index) = wid(index) + 3;
    [truefalse, aux] = ismember('FIT', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('dVFit', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('CFit', colheadings);wid(aux) = wid(index);
    [truefalse, aux] = ismember('Penal', colheadings);wid(aux) = wid(index);
    [truefalse, index] = ismember('N', colheadings);
    wid(index) = wid(index)+3;
    
    [truefalse, index] = ismember('lambda', colheadings);
    [truefalse, aux] = ismember('e', colheadings);wid(aux) = wid(index);

    wid(end) = wid(end) + 1;
    wid(end-1) = wid(end-1) + 1;
    wid(end-2) = wid(end-2) + 1;
    wid(end-3) = wid(end-3) + 1;
    
    dVFp = 100*delVFit/fitness;
    MvFp = 100*MvecFit/fitness;
    CFp = 100*CostFit/fitness;

    data1 = [FITNESS  fPenal delVFit MvecFit CostFit dVFp MvFp CFp N,...
        lambda e delV Mvec Cost];
    data2 = [pPL*100 Mpay mf mp mE];
    data = [data1 data2];

    displaytable(data,colheadings,wid,fms,rowheadings,fid);
    fprintf(fid,'\n');
    fclose('all');
end

