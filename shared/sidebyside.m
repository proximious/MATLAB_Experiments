% script to create side-by-side plots

figure; % initialize figure
subplot(1,2,1)
hold on

subplot(1,2,2)
hold on

x=0:0.1:10;

for n=[0,1,2,3]
    subplot(1,2,1) % say which plot you are editing
    y1=n*x;
    plot(x,y1)
    
    subplot(1,2,2) % now we are editing the second plot
    y2=n*x.^2;
    plot(x,y2)
end
subplot(1,2,1) % to change the ylimits on the first plot
ylim([0,50]) 
subplot(1,2,2) % to change the ylimits on the second plot
ylim([0,100])