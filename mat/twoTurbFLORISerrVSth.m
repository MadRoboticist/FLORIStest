close()
fontsize = 16;
x=(0:1:30);
l=10;
colors = {[0, 0.4470, 0.7410],...
          [1, 0, 0],...	          	
          [0.5, 1, 0.1250],	  ...        	
          [1, 0.1840, 1],	     ...     	
          [0, 0.75, 1],	        ...  	
          [0.75, 0.7450, 0.0],	          ...	
          [0.6350, 0.0780, 0.1840],...
          [0, 0, 0],...
          [0,0.5,0]}
SPDerr = subplot(2,1,1);
%title({'Wind Speed Estimation Error','Convergence by \gamma'});
set(get(SPDerr,'XLabel'), 'String', 'iterations');
set(get(SPDerr,'YLabel'), 'String', '$e_v$ (m/s)','Interpreter','latex');
xlim([0,l]);
hold on
plot(x, spdErrTh0,'color',colors{1})
plot(x, spdErrTh1,'--','color',colors{2})
plot(x, spdErrTh2,'color',colors{3})
plot(x, spdErrTh3,'color',colors{4})
plot(x, spdErrTh4,'color',colors{5})
plot(x, spdErrTh5,'color',colors{6})
plot(x, spdErrTh6,'--','color',colors{7})
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
plot(x, rad2deg(dirErrTh0),'color',colors{1})
plot(x, rad2deg(dirErrTh1),'--','color',colors{2})
plot(x, rad2deg(dirErrTh2),'color',colors{3})
plot(x, rad2deg(dirErrTh3),'color',colors{4})
plot(x, rad2deg(dirErrTh4),'color',colors{5})
plot(x, rad2deg(dirErrTh5),'color',colors{6})
plot(x, rad2deg(dirErrTh6),'--','color',colors{7})
%plot(x, rad2deg(dirErrTh7), '-.')
%plot(x, rad2deg(dirErrTh8), '-.')
%plot(x, rad2deg(dirErrTh9), '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
legend('no mask','\gamma=0.1','\gamma=0.2','\gamma=0.3','\gamma=0.4','\gamma=0.5','\gamma=0.6','\gamma=0.7','\gamma=0.8','\gamma=0.9')
