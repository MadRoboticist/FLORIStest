close()
fontsize = 14;
x=(0:1:200);
l=60;
SPDerr = subplot(2,1,1);
title({'Wind Speed Estimation Error','Convergence by \gamma'});
set(get(SPDerr,'XLabel'), 'String', 'iterations');
set(get(SPDerr,'YLabel'), 'String', 'error (m/s)');
xlim([0,l]);
hold on
plot(x, th0spd)
plot(x, th1spd)
plot(x, th2spd)
plot(x, th3spd)
plot(x, th4spd)
plot(x, th5spd)
plot(x, th6spd)
plot(x, th7spd, '-.')
plot(x, th8spd, '-.')
plot(x, th9spd, '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
DIRerr = subplot(2,1,2);
title({'Wind Direction Estimation Error','Convergence by \gamma'})
set(get(DIRerr,'XLabel'), 'String', 'iterations')
set(get(DIRerr,'YLabel'), 'String', 'error (degrees)')
xlim([0,l])
hold on
plot(x, rad2deg(th0dir))
plot(x, rad2deg(th1dir))
plot(x, rad2deg(th2dir))
plot(x, rad2deg(th3dir))
plot(x, rad2deg(th4dir))
plot(x, rad2deg(th5dir))
plot(x, rad2deg(th6dir))
plot(x, rad2deg(th7dir), '-.')
plot(x, rad2deg(th8dir), '-.')
plot(x, rad2deg(th9dir), '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
legend('no mask,','\gamma=0.1,','\gamma=0.2,','\gamma=0.3,','\gamma=0.4,','\gamma=0.5,','\gamma=0.6,','\gamma=0.7,','\gamma=0.8,','\gamma=0.9')
