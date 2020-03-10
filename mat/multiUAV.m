UAVs = {UAV0, UAV1, UAV2, UAV3, UAV4, UAV5, UAV6, UAV7, UAV8};
iters = [0:1:length(UAV0(:,1))-2];
V = zeros(1,length(UAV0(:,1))-1);
L=length(UAV0(:,1));
L=30
colors = {[0, 0.4470, 0.7410],...
          [1, 0, 0],...	          	
          [0.5, 1, 0.1250],	  ...        	
          [1, 0.1840, 1],	     ...     	
          [0, 0.75, 1],	        ...  	
          [0.75, 0.7450, 0.0],	          ...	
          [0.6350, 0.0780, 0.1840],...
          [0, 0, 0],...
          [0,0.5,0]}
for i=1:length(UAV0(:,1))-1
   for j=1:length(UAVs)
      V(i)=V(i)+UAVs{j}(i,1)^2+(UAVs{j}(i,2)*180/pi)^2;
   end
end
figure('DefaultAxesFontSize',16)
subplot(1,3,1)
for i=1:length(UAVs)
    plot(iters-1, UAVs{i}(1:L,1),'color',colors{i})
    hold on
end
title("Speed Error",'Interpreter','latex')
ylabel("m/s",'Interpreter','latex')
xlabel("iterations",'Interpreter','latex')
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5',...
    'Agent 6','Agent 7','Agent 8','Agent 9','Interpreter','latex')
grid on
subplot(1,3,2)
for i=1:length(UAVs)
    plot(iters-1, UAVs{i}(1:L,2)*180/pi,'color',colors{i})
    hold on
end
    
ylabel("degrees $(^\circ)$",'Interpreter','latex');
xlabel("iterations",'Interpreter','latex')
title("Direction Error",'Interpreter','latex')
grid on
subplot(1,3,3)
%semilogy(iters-1, V)
plot(iters-1,V)
ylabel("$\displaystyle\sum_{i=1}^9 V(\hat{x_i})$",'Interpreter','latex')
xlabel("iterations",'Interpreter','latex')
title("Lyapunov Function",'Interpreter','latex')
grid on