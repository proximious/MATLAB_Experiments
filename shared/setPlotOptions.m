function setPlotOptions()
% to reset back to factory defaults, enter into the Command Window: reset(groot) 

    set(groot, 'DefaultFigureColor', 'w'); % White background for figures
    set(groot, 'DefaultAxesFontSize', 12); % Font size for axes labels
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.5); % Scale for axes labels
    set(groot, 'DefaultAxesFontWeight', 'normal'); % Normal/Bold font for axes labels
    set(groot, 'DefaultLineLineWidth', 2); % Default line width
    set(groot, 'DefaultAxesLineWidth', 1); % Axes line width
    set(groot, 'DefaultAxesGridLineStyle', '--'); % Dashed grid lines
    set(groot, 'DefaultAxesBox', 'on'); % Box around plots
    set(groot, 'DefaultAxesTickLength', [0.02, 0.02]); % Tick sizes
    set(groot, 'DefaultAxesTickDir', 'in'); % Tick direction

    % Enable grid by default, if you want
    set(groot, 'DefaultAxesXGrid', 'off');
    set(groot, 'DefaultAxesYGrid', 'off');

    % Display a message indicating that plot options have been set
    disp('Global plot options set');
end