function val = fitn(x)
global cut Mpay W1 W2 e Udes DEBUG

W1 = 1;
W2 = 1000000;


e = 0.1;
out = divVec(x,cut);

N = de2re(out{1},1,5);

lambda = de2re(out{2},0.0001,1);

obj(1) =  abs(real(log(lambda + e*(1 - lambda)))^(-N)/W1);

obj(2) = (Mpay/((lambda)^N));

val = obj(1) + obj(2)/W2;

if val > 10; val = 10;end

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
