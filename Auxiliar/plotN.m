function plotN(x, flag)
global g0 Isp cut test type out finalName PATH
global hMarks

switch flag
    case 'init'
        e = [0.001 0.01 0.05 0.1 0.2 0.5];
        l = 0.0001:0.0001:1;
        
        figure('Name','plotN','NumberTitle','off')
        for i=1:size(e,2),
            U=g0*Isp*log((1+l)./(l+e(i)));
            dataPlotN = semilogx(l,U);
            hold on
        end
        hMarks=[];
        legend('e = 0.001','e = 0.01','e = 0.05','e = 0.1','e = 0.2','e = 0.5')
        
    case {'iter','interrupt'}
        
        if strcmp(type,'bitString'),
            
            out = divVec(x,cut);
            N = de2re(out{1},test.Nmin,test.Nmax);
            lBest = de2re(out{2},test.Lmin,test.Lmax);
            eBest = de2re(out{3},test.Emin,test.Emax);
            
        elseif strcmp(type,'doubleVector'),
            N = round(x(1));
            eBest = x(3);
            lBest = x(2);
        end
        U=g0*Isp*log((1+lBest)./(lBest+eBest));
        
        f = findobj('name','plotN');
        set(0, 'currentfigure', f);
        hold on
        set(hMarks,'Color','r','LineWidth',1, 'Marker','o','MarkerSize',5);
        
        hMarks = [hMarks semilogx(lBest,U,'kx', 'LineWidth',3,'MarkerSize',9)];
        
        drawnow
        
%         newX = [get(dataPlotN,'Xdata') lBest];
%         newY = [get(dataPlotN,'Ydata') U];
%         set(dataPlotN,'Xdata',newX, 'Ydata',newY);
%         semilogx(lBest,U,'bo')

        hold off
        
        
    case 'done'
        
        saveas(gcf, strcat(fullfile(PATH, finalName),'.fig'))
end
end
