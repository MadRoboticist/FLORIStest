UAVs = {UAV0, UAV1, UAV2, UAV3, UAV4, UAV5, UAV6, UAV7, UAV8};
iters = [1:1:length(UAV0(:,1))-1];
V = zeros(1,length(UAV0(:,1))-1);

for i=1:length(UAV0(:,1))-1
   for j=1:length(UAVs)
      V(i)=V(i)+UAVs{j}(i,1)^2+(UAVs{j}(i,2)*180/pi)^2;
   end
end
subplot(1,3,1)
for i=1:length(UAVs)
    plot(iters, UAVs{i}(1:30,1))
    hold on
end
subplot(1,3,2)
for i=1:length(UAVs)
    plot(iters, UAVs{i}(1:30,2)*180/pi)
    hold on
end
subplot(1,3,3)
plot(iters, V)