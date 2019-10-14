close()
fontsize = 16;
x=(0:1:30);
l=10;
SPDerr = subplot(2,1,1);
%title({'Wind Speed Estimation Error','Convergence by \gamma'});
set(get(SPDerr,'XLabel'), 'String', 'iterations');
set(get(SPDerr,'YLabel'), 'String', '$e_v$ (m/s)','Interpreter','latex');
xlim([0,l]);
hold on
plot(x, spdErrTh0)
plot(x, spdErrTh1)
plot(x, spdErrTh2)
plot(x, spdErrTh3)
plot(x, spdErrTh4)
plot(x, spdErrTh5)
plot(x, spdErrTh6)
%plot(x, spdErrTh7, '-.')
%plot(x, spdErrTh8, '-.')
%plot(x, spdErrTh9, '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
DIRerr = subplot(2,1,2);
%title({'Wind Direction Estimation Error','Convergence by \gamma'})
set(get(DIRerr,'XLabel'), 'String', 'iterations')
set(get(DIRerr,'YLabel'), 'String', '$e_\theta$ ($^\circ$)','Interpreter','latex')
set(gca,'fontsize', fontsize);
xlim([0,l])
hold on
plot(x, rad2deg(dirErrTh0))
plot(x, rad2deg(dirErrTh1))
plot(x, rad2deg(dirErrTh2))
plot(x, rad2deg(dirErrTh3))
plot(x, rad2deg(dirErrTh4))
plot(x, rad2deg(dirErrTh5))
plot(x, rad2deg(dirErrTh6))
%plot(x, rad2deg(dirErrTh7), '-.')
%plot(x, rad2deg(dirErrTh8), '-.')
%plot(x, rad2deg(dirErrTh9), '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
legend('no mask','\gamma=0.1','\gamma=0.2','\gamma=0.3','\gamma=0.4','\gamma=0.5','\gamma=0.6','\gamma=0.7','\gamma=0.8','\gamma=0.9')
