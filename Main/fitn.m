function obj = fitn(x)
global cut Mpay

%cantileverVolume Calculate volume of a stepped cantilever
%
%   V = cantileverVolume(x) calculates the volume of the stepped cantilever
%   in the "Solving a Mixed Integer Engineering Design Problem Using the
%   Genetic Algorithm" example. 

%   Copyright 2012 The MathWorks, Inc.

% Volume of cantilever beam

g=sprintf('%d ', x);
fprintf('Answer: %s\n', g)
e = 0.1;
out = divVec(x,cut);

N = de2re(out{1},1,5);
fprintf('N: %f \n',N)
lambda = de2re(out{2},0.0001,1);
fprintf('lambda: %f \n',lambda)
obj = abs(real(log(lambda + e*(1 - lambda)))^(-N));
%obj(1) =  abs(real(log(lambda + e*(1 - lambda)))^(-N));
%obj(2) = Mpay/((lambda)^N)
%U = 100*(x(1)*x(2) + x(3)*x(4) + x(5)*x(6) + x(7)*x(8) + x(9)*x(10));
