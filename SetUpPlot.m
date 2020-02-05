function [fig, Splot, Hplot, Tplot, Xplot, HXoft, plotcount] = SetUpPlot(S,H,X,C,xname,x2plot,HXoft,plotcount,K,J,dx,dy,relief,vexag)

% NaN out the channels to make plots look better.
H(C==1) = NaN;
tempX1 = X(:,:,1);
tempX2 = X(:,:,2);
tempX1(C==1) = NaN;
tempX2(C==1) = NaN;
X(:,:,1) = tempX1;
X(:,:,2) = tempX2;

fig=figure; figure(fig)
set(fig, 'color', 'w');
colormap jet

Splot = subplot(2,2,1);
        vertexag = [1 1 1]; 
        vertexag(2) = K/J; % ratio of x to y axis length, next line is vertical exaggeration
        vertexag(3) = relief*vexag/(J*dx);
        axes(Splot)
        surf(S); shading flat
        title('Surface elevation (m). Iteration=0, t=0')
        shading flat % other options are flat, interp
        axis([0 J 0 K 0 relief])
        set(Splot, 'PlotBoxAspectRatioMode', 'Manual', 'PlotBoxAspectRatio', vertexag, 'nextplot', 'replacechildren');
        
Hplot = subplot(2,2,2);
        axes(Hplot)
        surf(H); shading flat
        title('Soil thickness (m)')
        set(Hplot, 'nextplot', 'replacechildren');

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

Xplot = subplot(2,2,4);
        axes(Xplot)
        surf(X(:,:,x2plot)); shading flat
        title(['Mineral concentration (' num2str(x2plot) ', ' xname{x2plot} ')'])
        set(Xplot, 'nextplot', 'replacechildren');

drawnow

plotcount = plotcount+1;