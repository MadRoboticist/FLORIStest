close()
fontsize=14;
per_covg_max = [0.0,2,5,10:10:100]
gamma = (0:10:100)
iters = [NaN,5, 5, 5, 5, 5, 5,...
         5, 5, 4, 4, 4, 4]
mask_max = [ 0.0, 1.8, 4.86, 9.86, 17.77, 24.58, 31.38,...
     43.33, 45.41, 47.36,53.61,56.67,56.94];

mask_min = [ 0.0, 1.11, 4.17, 8.75, 14.72, 22.36, 26.11,...
             32.08, 36.81, 41.67,46.25,49.86,48.75];
 
mask_avg = [ 0.0, 1.41, 4.54, 9.26, 16.02, 23.33, 28.47,...
             40.07, 43.1, 44.69,51.44,54.19,53.25];
         
max_per = [100,31.50, 26.25, 20.97, 18.10,14.02,8.89,7.08,3.47,2.22,0];
min_per = [100,25.83,20.69,14.44,10.97,6.25,2.08,0.97,0.28,0.14,0];
avg_per = [100,28.80,23.52,18.19,14.13,9.82,6.52,3.70,1.74,0.32,0];
converg = [5,9,11,14,17,23,34,61,135,400,100000000];       
yyaxis right
plot(per_covg_max,iters,'-.o');

ylim([0,10])
ylabel('steps to convergence');
hold on
%plot(fliplr(gamma),converg,'-.o')
yyaxis left
plot(per_covg_max,mask_max,'-*',per_covg_max,mask_min,'-x',per_covg_max,mask_avg,'-o')
ylabel('% coverage');
xlabel('$\Gamma$','Interpreter','latex');
legend('max. coverage','min. coverage','avg. coverage','convergence')
ylim([0, 100])
hold on
%plot(fliplr(gamma),max_per,'--*',fliplr(gamma),min_per,'--x',fliplr(gamma),avg_per,'-o')
set(gca,'fontsize', fontsize);
set(gca, 'XDir','reverse')