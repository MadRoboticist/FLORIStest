close()
fontsize = 16;
l=10;
x=(0:1:10);
SPDerr = subplot(2,1,1);
%title({'Wind Speed Estimation Error','Convergence by M'});
set(get(SPDerr,'XLabel'), 'String', 'iterations');
set(get(SPDerr,'YLabel'), 'String', '$e_v$ (m/s)','Interpreter','latex');
xlim([0,l]);
hold on
plot(x, UAV02(:,1))
plot(x, UAV05(:,1))
plot(x, UAV1(:,1))
plot(x, UAV2(:,1))
plot(x, UAV3(:,1))
plot(x, UAV4(:,1))
plot(x, UAV5(:,1))
plot(x, UAV6(:,1))
plot(x, UAV7(:,1), '-.')
plot(x, UAV8(:,1), '-.')
plot(x, UAV9(:,1), '-.')
plot(x, UAV10(:,1), '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
DIRerr = subplot(2,1,2);
%title({'Wind Direction Estimation Error','Convergence by M'})
set(get(DIRerr,'XLabel'), 'String', 'iterations')
set(get(DIRerr,'YLabel'), 'String', '$e_\theta$ ($^\circ$)','Interpreter','latex')
xlim([0,l])
hold on
plot(x, UAV02(:,2))
plot(x, UAV05(:,2))
plot(x, UAV1(:,2))
plot(x, UAV2(:,2))
plot(x, UAV3(:,2))
plot(x, UAV4(:,2))
plot(x, UAV5(:,2))
plot(x, UAV6(:,2))
plot(x, UAV7(:,2), '-.')
plot(x, UAV8(:,2), '-.')
plot(x, UAV9(:,2), '-.')
plot(x, UAV10(:,2), '-.')
grid minor
grid on
hold off
set(gca,'fontsize', fontsize);
legend('$\Gamma=2\%$','$\Gamma=5\%$','$\Gamma=10\%$','$\Gamma=20\%$','$\Gamma=30\%$','$\Gamma=40\%$',...
    '$\Gamma=50\%$','$\Gamma=60\%$','$\Gamma=70\%$','$\Gamma=80\%$','$\Gamma=90\%$','$\Gamma=100\%$','Interpreter','latex')
