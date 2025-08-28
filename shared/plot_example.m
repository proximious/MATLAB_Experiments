close all

x=-10:0.01:10;
y=sin(x);


plot(x,y)
xlabel('$x$','interpreter','latex')
ylabel('$\sin x$','interpreter','latex')

xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
xticklabels({'-3\pi','-20\pi','-\pi','0','\pi','2\pi','big number'})
% note you can mislabel the points on purpose!

yticks([-1 -0.8 -0.2 0 0.2 0.8 1])