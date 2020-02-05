function [HXoft, plotcount] = DrawPlot(fig,Splot,Hplot,Tplot,Xplot,S,H,X,C,xname,x2plot,HXoft,plotcount,n,t)

% NaN out the channels to make plots look better.
H(C==1) = NaN;
tempX1 = X(:,:,1);
tempX2 = X(:,:,2);
tempX1(C==1) = NaN;
tempX2(C==1) = NaN;
X(:,:,1) = tempX1;
X(:,:,2) = tempX2;

figure(fig)

axes(Splot)
        surf(S);         
        shading flat
        title(['Surface elevation (m). Iteration=' num2str(n) ', t=' num2str(t)])

axes(Hplot)
        surf(H); shading flat
        title('Soil thickness (m)')
        
Tplot = subplot(2,2,3);
        HXoft(plotcount,2) = nanmean(H(:)); 
        HXoft(plotcount,3) = nanmean(nanmean(X(:,:,x2plot)));  
        [Tplot hH hX] = plotyy(HXoft(1:plotcount,1),HXoft(1:plotcount,2),HXoft(1:plotcount,1),HXoft(1:plotcount,3));
        
        title(['H(t) and X(t) (' num2str(x2plot) ', ' xname{x2plot} ')'])
        xlabel('Time (yr)') 
        
        set(get(Tplot(1),'Ylabel'),'String','Mean soil thickness (m)') 
        set(get(Tplot(2),'Ylabel'),'String',['Mean concentration (' num2str(x2plot) ', ' xname{x2plot} ')']) 
        set(Tplot(1),'xlim',[0 max(HXoft(:,1))]);        
        set(Tplot(2),'xlim',[0 max(HXoft(:,1))]);        
        
axes(Xplot)
        surf(X(:,:,x2plot)); shading flat
        title(['Mineral concentration (' num2str(x2plot) ', ' xname{x2plot} ')'])

drawnow

plotcount = plotcount + 1;