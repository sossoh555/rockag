function plotN
e = [0.001 0.01 0.05 0.1 0.2 0.5];
l = 0.0001:0.0001:1;


for i=1:size(e,2),
U=log((1+l)./(l+e(i)));
    semilogx(l,U)
    hold on
end

end
