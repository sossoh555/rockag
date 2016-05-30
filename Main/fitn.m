function FITNESS = fitn(x)
global cut Mpay g0 Isp Udes DEBUG type test POP
ALLbool = false;
MVECbool = true;
DELVbool = true;
COSTbool = true;


A = 0.6:-0.1:0.1;
B = 15:2:25;
W1 = 1;
W2 = 1E6;
W3 = 4E2;
Cost = 0;delV = 0;Mvec = 0; mf = 0; mp = 0; mE = 0; 
CostFit = 0;delVFit = 0;MvecFit = 0;

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
%Important Considerations
%=======================
pPL = lambda/(1+lambda);
    
    
%=======================
%          U
%=======================

if DELVbool || ALLbool,
   % U = log(((1+lambda)/(e + lambda))^N);
   %delV = g0*Isp*U;
    
    num = (pPL^(1/N))*(1-e) + e;
    U = log(1/num);
    delV = U*Isp*g0*N;
    
    delVFit = abs((Udes - delV)/Udes/W1);
end




%=======================
%         Mvec
%=======================

if MVECbool || ALLbool,
    mf = Mpay*((lambda+e)/lambda);
    num = (pPL^(1/N))*(1-e) + e;
    n = (1/num)^N;
    Mvec = mf*n;
    MvecFit = Mvec/W2;
end

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
if COSTbool || ALLbool,
    Nr = round(N);
    Cost = 0;
    for i=1:Nr,
        Cost = Cost + Nr*(1 + A(Nr)*exp(-((0.025*B(Nr))/(lambda-lambda*e+e))^3.5));
    end
    CostFit = Cost/W3;
end


%=======================
%        FITNESS
%=======================
FITNESS = delVFit + MvecFit + CostFit;
if FITNESS >1, FITNESS = 1; end

if DEBUG == false
%nameFile =  strcat(POP,finalName,'.txt');
%fid = fopen(strcat(fullfile(PATH, nameFile)), 'at' );
fid = fopen(POP, 'at' );

            fprintf(fid,'%g ',x(:));
            fprintf(fid,'| %f',FITNESS);
            fprintf(fid,'\n');
    %%%%%%%%%%%%%%%% INFORMACOES DO AG %%%%%%%%%%%%%%%%%%%%%%%
    mp = Mvec - mf;
    mE = Mvec - mp - Mpay;
 
    ch0 = {'AllB','MvB','dVB','CB'};
    ch1 = {'FIT','dVFit','MvecFit','CFit','dVF [%]','MvF [%]','CF [%]','N','lambda','e','delV [kg/s]','Mvec [kg]','Cost [U$]'};
    ch2 = {'pPL [%]','Mpay [kg]','mf [kg]','mp [kg]','mE [kg]'};
    
    colheadings = [ch0 ch1 ch2];
    
    fms0 = {'d','d','d','d'};
    fms1 = {'.2E','.2E','.2E','.2E','.3f','.3f','.3f','.2f','.4f','.4f','.4f','.2E','.4f'};
    fms2 = {'.4f','.2E','.2E','.2E','.2E'};
    
    fms = [fms0 fms1 fms2];
    
    rowheadings = {''};
    
    wid = ones(1,size(colheadings,2));
    for i=1:size(colheadings,2),
        wid(i) = length(colheadings{1,i}) ; 
    end
    wid(7) = wid(7) + 1;
    wid(5) = wid(7); wid(6) = wid(7); wid(8) = wid(7); wid(14) = wid(13);wid(12) = wid(12)+3;
    
    wid(end) = wid(end) + 1;
    wid(end-1) = wid(end-1) + 1;
    wid(end-2) = wid(end-2) + 1;
    wid(end-3) = wid(end-3) + 1;
    dVFp = 100*delVFit/FITNESS;
    MvFp = 100*MvecFit/FITNESS;
    CFp = 100*CostFit/FITNESS;
    
    data0 = [ALLbool MVECbool DELVbool COSTbool];
    data1 = [FITNESS delVFit MvecFit CostFit dVFp MvFp CFp N lambda e delV Mvec Cost];
    data2 = [pPL*100 Mpay mf mp mE];
    data = [data0 data1 data2];
    

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

