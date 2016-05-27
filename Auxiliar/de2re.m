
function val = de2re(x, lb, ub)
global DEBUG
re =  bi2de(x, 'left-msb');
val = re*(ub - lb)/(2^length(x) - 1) + lb;

if DEBUG fprintf('bi2de: %f \t value: %f \n',re,val)
    
end
