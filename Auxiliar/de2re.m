
function val = de2re(x, lb, ub)
global DEBUG
re =  bi2de(x, 'left-msb');
val = re*(ub - lb)/(2^length(x) - 1) + lb;

if DEBUG fprintf('\nbi2de: %.0f \t value: %.4f',re,val)
    
end
