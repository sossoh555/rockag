function val = fitn(x)
global cut Mpay W1 W2 g0 Isp Udes DEBUG

A = 0.6:-0.1:0.1;
B = 15:2:25;

W2 = 200000;
lambda = x(1);
e = x(2);
%N = round(x(3));

%=======================
%          U
%=======================
U = log((1+lambda)/(e + lambda));
delV = abs(Udes-g0*Isp*U)/Udes;


%=======================
%         Mvec
%=======================
pPL = lambda/(1+ lambda);


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
%Cost = 


%=======================
%         Val
%=======================
Mvecfit =0;
val = delV; %+ Mvecfit;
fprintf('\nLamb  \t estr \t DelvaV \t Mvec \t\t Cost')
fprintf('\n%f  \t %f \t %f \t %f \t %f\n\n',lambda,e,delV)
% W1 = 1;
% W2 = 1000000;
% 
% 
% 
% %out = divVec(x,cut);
% 
% %N = de2re(out{1},1,5);
% N = x(1);
% %lambda = de2re(out{2},0.0001,1);
% lambda = x(2);
% 
% e = 0.1;
% 
% 
% obj(1) =  abs(Udes - real(log(lambda + e*(1 - lambda))^(N)));
% 
% obj(2) = (Mpay/((lambda)^N));
% 
% val = obj(1)/W1 + obj(2)/W2;
% 
% if val > 5; val = 5;end

if DEBUG
g = sprintf('%d ', x);
fprintf('Answer: %s\n', g)
fprintf('N: %f \n',N)
fprintf('lambda: %f \n',lambda)
fprintf('obj(1): %f \n',obj(1))
fprintf('obj(2): %f \n',obj(2))
fprintf('val: %f \n',val)
fprintf('\n\n')
end

%fprintf('\n\n',obj(2))

%U = 100*(x(1)*x(2) + x(3)*x(4) + x(5)*x(6) + x(7)*x(8) + x(9)*x(10));
