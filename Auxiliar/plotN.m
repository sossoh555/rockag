function plotN(x, flag)
global g0 Isp cut test type out finalName PATH
global hMarks Udes

fprintf('flag: %s\n',flag)
e = [0.001 0.01 0.05 0.1 0.2 0.5];
l = 0.0001:0.0001:1;
switch flag
    case 'init'
        
        N = 1;
        
        figure('Name','plotN','NumberTitle','off')
        for i=1:size(e,2),
            [pPL, U] = delVcalc(l,e(i),N);
            %U=g0*Isp*log((1+l)./(l+e(i)));
            semilogx(l,U,'Tag',sprintf('plotN%d',i));
            hold on
        end

        
        hMarks=[];
        title('\Delta V vs \lambda para N = 1');
        xlabel('\lambda');
        ylabel('\Delta V [km/s]');
        
        
    case {'iter','interrupt'}
        
        if strcmp(type,'bitString'),
            
            out = divVec(x,cut);
            N = de2re(out{1},test.Nmin,test.Nmax);
            lBest = de2re(out{2},test.Lmin,test.Lmax);
            eBest = de2re(out{3},test.Emin,test.Emax);
            
        elseif strcmp(type,'doubleVector'),
            N = x(1);
            eBest = x(3);
            lBest = x(2);
        end
        
        [pPL, U] = delVcalc(lBest,eBest,N);
        %U=g0*Isp*log((1+lBest)./(lBest+eBest));
        
        f = findobj('name','plotN');
        set(0, 'currentfigure', f);
        hold on
        
        
        axesO = get(gca,'children');
        if isempty(findobj('Marker','x'))
            semilogx(lBest,U,'kx', 'LineWidth',3,'MarkerSize',9);
        else
            if isempty(findobj('Marker','o'))
                plotBest = findobj(get(gca,'Children'),'Marker','x');
                newX = [get(plotBest,'Xdata')];
                newY = [get(plotBest,'Ydata')];
                semilogx(newX,newY,'ro', 'LineWidth',1,'MarkerSize',5);
                chH = get(gca,'Children');
                if chH(2).Marker == 'x'
                    set(gca,'Children',[chH(2);chH(1); chH(3:end)])
                    legend('e = 0.001','e = 0.01','e = 0.05','e = 0.1','e = 0.2','e = 0.5',...
                        'Melhores indivíduos anteriores','Melhor indivíduo atual');
                end
            else
                          plotBest = findobj(get(gca,'Children'),'Marker','x');
                plotOld = findobj(get(gca,'Children'),'Marker','o');

                newX = [get(plotOld,'Xdata') get(plotBest,'Xdata')];
                newY = [get(plotOld, 'Ydata') get(plotBest,'Ydata')];
                set(plotOld, 'Xdata',newX, 'Ydata', newY);
                
       
                
                
              
            end
            plotBest = findobj(get(gca,'Children'),'Marker','x');
            newX = [lBest];
            newY = [U];
            set(plotBest, 'Xdata',newX, 'Ydata', newY);
            %set(hMarks,'Color','r','LineWidth',1, 'Marker','o','MarkerSize',5);
            
            
            
            %hMarks = [hMarks semilogx(lBest,U,'kx', 'LineWidth',3,'MarkerSize',9)];
            
        end
        
        
        drawnow
        
        
        %         newX = [get(dataPlotN,'Xdata') lBest];
        %         newY = [get(dataPlotN,'Ydata') U];
        %         set(dataPlotN,'Xdata',newX, 'Ydata',newY);
        %         semilogx(lBest,U,'bo')
        
    case 'done'
        
             
%              out = divVec(x,cut);
%             N = de2re(out{1},test.Nmin,test.Nmax);
%             lBest = de2re(out{2},test.Lmin,test.Lmax);
%             eBest = de2re(out{3},test.Emin,test.Emax);
%        f = findobj('name','plotN');
%         set(0, 'currentfigure', f);
%         hold on
%                   for i=1:size(e,2),
%             [pPL, U] = delVcalc(l,e(i),6);
%             %U=g0*Isp*log((1+l)./(l+e(i)));
%             semilogx(l,U);
%             hold on
%         end
%         
                
        % f = findobj('name','plotN');
        % set(0, 'currentfigure', f);
        % legend('e = 0.001','e = 0.01','e = 0.05','e = 0.1','e = 0.2','e = 0.5','Previous Best','Current Best');
        % drawnow
        
        %saveas(gcf, strcat(fullfile(PATH, finalName),'.fig'))
end
end
