close all

x=linspace(-3,3,50);
y1=x.^2-1;
y2=1-x.^2;


plot(x,y1,'r-')
hold on
plot(x,y2,'b--')

% the annotation command 
% Annoyingly, these are not given in actual coordinates.
% Instead, it is the fraction of the x domain, 0->1, and similar for y.
% So you can tweak until it looks good, which is what I usually do,
% or try to be more systematic about it by, e.g., setting xlim as
% xlim(min(x),max(x))
% then the position of the arrow tip, in a fraction of the domain, is
% (x_intersect-min(x))/(max(x)-min(x))
%
% (Honestly I usually just tweak it and re-run until it looks good.)

annotation('textarrow',[0.28,0.38],[0.518,0.518],'String','Intersect 1')
annotation('textarrow',[0.76,0.66],[0.518,0.518],'String','Intersect 2')

% [0.2333,0.39],[0.52,0.52] are [x_start,x_end] and [y_start,y_end]
% I am not sure why y_end ~= 0.5 doesn't look right hmm