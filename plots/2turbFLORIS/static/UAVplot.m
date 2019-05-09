close()
fontsize = 14;
iters=(0:1:150)
dir = subplot('Position',[0.09 0.1 0.75 0.35])
spd = subplot('Position',[0.09 0.55 0.75 0.35])
for i=1:1:7
    subplot(spd)
    plot(iters,evalin('base',sprintf('spdErrUAV%d',i*5)),'-.')
    
    hold on
    subplot(dir)
    plot(iters,rad2deg(evalin('base',sprintf('dirErrUAV%d',i*5))),'-.')
    hold on
end
for i=8:1:14
    subplot(spd)
    plot(iters,evalin('base',sprintf('spdErrUAV%d',i*5)),'--')
    
    hold on
    subplot(dir)
    plot(iters,rad2deg(evalin('base',sprintf('dirErrUAV%d',i*5))),'--')
    hold on
end
for i=15:1:20
    subplot(spd)
    plot(iters,evalin('base',sprintf('spdErrUAV%d',i*5)))
    
    hold on
    subplot(dir)
    plot(iters,rad2deg(evalin('base',sprintf('dirErrUAV%d',i*5))))
    hold on
end
subplot(spd); xlim([0,60]); title({'Speed Estimation Error Convergence','by Max. Percent Coverage'});
set(gca,'fontsize', fontsize); grid minor; grid on; ylim([0,6]); xlim([0,61]);
ylabel('m/s');
subplot(dir); xlim([0,60]); title({'Direction Estimation Error Convergence','by Max. Percent Coverage'});
set(gca,'fontsize', fontsize); grid minor; grid on; ylim([0,12]); xlim([0,61]);
ylabel('degrees'); xlabel('iterations');
for i=1:20
    titl{i}= num2str(i*5);
    titl{i}=strcat(titl{i},'%');
end
linkaxes([spd,dir],'x')
legend(titl,'Interpreter','tex','Position',[0.85 0.05 0.15 0.9])