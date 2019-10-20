close all; clc;
iters = [1:1:length(ERRnorm0)];
mean(ERRnorm0)
mean(refERRnorm0)
figure(); fontsize = 16;
plot(iters, ERRnorm0,'r');
hold on;
plot(iters, refERRnorm0,'b');
plot([0 length(ERRnorm0)],[mean(ERRnorm0) mean(ERRnorm0)],'--r');
plot([0 length(ERRnorm0)],[mean(refERRnorm0) mean(refERRnorm0)],'--b');
xlabel("iterations");
ylabel("$||e_{u-field}||$",'Interpreter','latex');
xlim([0 length(ERRnorm0)])
legend("estimated $u-$field","static $u-$field",...
        "mean estimate error","mean static error",...
        'Interpreter','latex');
title({"Norms of error fields between FLORIS and SOWFA",...
        "static reference vs. moving estimate"});
set(gca,'FontSize',fontsize);