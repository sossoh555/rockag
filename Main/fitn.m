function val = fitn(x)
global cut Mpay W1 W2 e

W1 = 1;
W2 = 100000;

g = sprintf('%d ', x);
fprintf('Answer: %s\n', g)
e = 0.1;
out = divVec(x,cut);

N = de2re(out{1},1,5);
fprintf('N: %f \n',N)

lambda = de2re(out{2},0.0001,1);
fprintf('lambda: %f \n',lambda)

obj(1) =  abs(real(log(lambda + e*(1 - lambda)))^(-N)/W1);
fprintf('obj(1): %f \n',obj(1))

obj(2) = (Mpay/((lambda)^N));
fprintf('obj(2): %f \n',obj(2))

val = obj(1) + obj(2)/W2;


%U = 100*(x(1)*x(2) + x(3)*x(4) + x(5)*x(6) + x(7)*x(8) + x(9)*x(10));
