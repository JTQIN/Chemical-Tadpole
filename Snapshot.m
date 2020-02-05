function Snapshot(S,n,t,run_name)

finalfig = figure;

% map view with curvature as color - requires curvn.m to calculate
% normalized curvature

imagesc(curvn(S)); axis image

title(['Iteration=' num2str(n) ', t=' num2str(t)])
drawnow              

print(finalfig, '-dpng', run_name)
close(finalfig)
