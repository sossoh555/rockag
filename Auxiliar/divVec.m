%http://stackoverflow.com/questions/29860008/split-vector-in-matlab

function arr  = divVec(vec,cut)
k=1;j=1;
for i=1:1:length(vec)
    if cut(i)==0
        arr{k}(j)=vec(i);
        j=j+1;
    else
        arr{k}(j)=vec(i);
        k=k+1;
        j=1;
    end
end
end