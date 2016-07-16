function plotN(x, flag)
global cut test type out
global hMarks Udes

fprintf('flag: %s\n',flag)
e = [0.001 0.01 0.05 0.1 0.2 0.5];
l = 0.00001:0.0001:1;
switch flag
    case 'init'
        
        N = 1;
        figure('Name','plotN','NumberTitle','off')
        
        for i=1:size(e,2),
            [pPL, U] = delVcalc(l,e(i),N);
            
            semilogx(l,U,'Tag',sprintf('plotN%d',i));
            hold on
        end
        
        semilogx([l(1),l(end)],[Udes Udes], 'k-',...
            'LineWidth',1.0, 'Tag','Udes');
        hMarks=[];
        title('\Delta V vs \lambda para N = 1','interp','none');
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
        
        f = findobj('name','plotN');
        set(0, 'currentfigure', f);
        hold on
        
        for i=1:size(e,2),
            plotNdata = findobj(get(gca,'Children'),...
                'Tag',sprintf('plotN%d',i));
            delete(plotNdata)
        end
        
        for i=1:size(e,2),
            [pPL, Ug] = delVcalc(l,e(i),N);
            semilogx(l,Ug,'Tag',sprintf('plotN%d',i));
        end
        
        set(f,'defaulttextinterpreter','latex');
        title(sprintf('Delta V vs lambda para N = %.2f',N))
        
        
        if isempty(findobj('Marker','x'))
            semilogx(lBest,U,'kx', 'LineWidth',3,'MarkerSize',9,...
                'Tag','CurrBest');
        else
            if isempty(findobj('Marker','o'))
                plotBest = findobj(get(gca,'Children'),'Marker','x');
                newX = [get(plotBest,'Xdata')];
                newY = [get(plotBest,'Ydata')];
                semilogx(newX,newY,'ro', 'LineWidth',1,'MarkerSize',5);
            else
                plotBest = findobj(get(gca,'Children'),'Marker','x');
                plotOld = findobj(get(gca,'Children'),'Marker','o');
                
                newX = [get(plotOld,'Xdata') get(plotBest,'Xdata')];
                newY = [get(plotOld, 'Ydata') get(plotBest,'Ydata')];
                set(plotOld, 'Xdata',newX, 'Ydata', newY);

                Utex = sprintf('U_{desejado} = %.2f km/s',Udes);

                plotBest = findobj(get(gca,'Children'),'Marker','x');
                plotOld = findobj(get(gca,'Children'),'Marker','o');
                plotUdes = findobj(get(gca,'Children'),'Tag','Udes');
                plotN1 = findobj(get(gca,'Children'),'Tag','plotN1');
                plotN2 = findobj(get(gca,'Children'),'Tag','plotN2');
                plotN3 = findobj(get(gca,'Children'),'Tag','plotN3');
                plotN4 = findobj(get(gca,'Children'),'Tag','plotN4');
                plotN5 = findobj(get(gca,'Children'),'Tag','plotN5');
                plotN6 = findobj(get(gca,'Children'),'Tag','plotN6');
                
                chH = get(gca,'Children');
                if strcmp(chH(end).Marker,'x'),
                    set(gca,'Children',...
                        [chH(end-2);chH(end-1);chH(end);chH(1:end-3)])
                     
                    legend([plotN1,plotN2,plotN3,plotN4,plotN5,plotN6,...
                        plotUdes,plotOld,plotBest,...
                    ],...
                    'e = 0.001','e = 0.01','e = 0.05','e = 0.1',...
                    'e = 0.2','e = 0.5',...
                    Utex,'Melhores individuos anteriores',...
                    'Melhor individuo atual')
                elseif strcmp(chH(end-1).Marker,'x'),
                    set(gca,'Children',...
                        [chH(end-1);chH(end-2);chH(end);chH(1:end-3)])
                    
                    legend([plotN1,plotN2,plotN3,plotN4,plotN5,plotN6,...
                        plotOld,plotUdes,plotBest],...
                    'e = 0.001','e = 0.01','e = 0.05','e = 0.1',...
                    'e = 0.2','e = 0.5',...
                    Utex,'1Melhores individuos anteriores',...
                    'Melhor individuo atual')
                end
                legend([...
                    plotN1,plotN2,plotN3,plotN4,plotN5,plotN6,...
                    plotUdes,plotOld,plotBest],...
                    'e = 0.001','e = 0.01','e = 0.05','e = 0.1',...
                    'e = 0.2','e = 0.5',...
                    Utex,'Melhores individuos anteriores',...
                    'Melhor individuo atual')
            end
            plotBest = findobj(get(gca,'Children'),'Marker','x');
            newX = [lBest];
            newY = [U];
            set(plotBest, 'Xdata',newX, 'Ydata', newY);
            
        end

        drawnow    
    case 'done'
        
end

end
