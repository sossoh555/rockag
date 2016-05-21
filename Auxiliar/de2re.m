
function val = de2re(x, lb, ub)
global DEBUG
re =  bi2de(x);
val = re*(ub - lb)/(2^length(x) - 1) + lb;

if DEBUG fprintf('re: %f \n',re)
    
end
