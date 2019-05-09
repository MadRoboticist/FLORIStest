close()
fontsize=14
lambda = (0.0:0.01:1.0);
length(lambda)
itersL = [30,32,32,32,32,32,32,32,32,32,...
         32,32,32,32,32,32,32,32,32,32,...
         32,32,32,32,32,32,32,32,32,32,...
         32,32,32,32,32,32,32,32,32,32,...
         32,32,32,32,32,32,32,31,39,39,...
         39,39,39,39,34,34,34,34,31,35,...
         33,27,29,29,29,30,31,31,32,32,...
         31,31,27,27,32,32,30,41,31,29,...
         34,37,32,30,33,32,34,36,31,34,...
         33,34,40,39,43,44,50,61,79,118,10000000];
     itersLavg(1)=30;
for i=1:20
    itersLavg(i+1)=mean(itersL(i:i+4));
end
coverage = (0:5:100);
itersCVG = [1000000,103,49,32,22,16,13,11,10,...
            10,9,8,8,8,8,7,7,7,7,7,7];
figure();
subplot('Position',[0.2 0.2 0.6 0.6]);
line(coverage,itersCVG,'color','k')
ylim([0,120])
set(gca,'fontsize', fontsize);
ax1 = gca;
set(ax1,'XColor','k','YColor','k')


xlabel('percent coverage')
ylabel('iterations to convergence')

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','b','YColor','k');
line(lambda,itersL,'Parent',ax2,'color','b')
ylim([0,120])
xlimits = get(ax1,'XLim');
ylimits = get(ax1,'YLim');
xinc = (xlimits(2)-xlimits(1))/10;
yinc = (ylimits(2)-ylimits(1));
set(ax1,'XTick',[xlimits(1):xinc:xlimits(2)],...
        'YTick',[ylimits(1):yinc:ylimits(2)])
title('Convergence vs. \lambda at 15% max. coverage')
xlabel('\lambda')

set(gca,'fontsize', fontsize);
     
         
