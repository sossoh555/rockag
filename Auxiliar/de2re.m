function val = de2re(x, lb, ub)

re =  bi2de(x);
fprintf('re: %f \n',re)
val = re*(ub - lb)/(2^length(x) - 1) + lb;

end
