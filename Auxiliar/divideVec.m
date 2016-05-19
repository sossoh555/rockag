function arr  = divideVec(vec,cut)
k=1;j=1;
for i=1:1:length(vec)
    if cut(i)==0
        arr{n}(j)=vec(i);
        j=j+1;
    else
        arr{n}(j)=vec(i);
        n=n+1;
        j=1;
    end
end
end