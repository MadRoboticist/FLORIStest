close()
fontsize=14;
per_covg_max = [0.0,2,5,10:10:100]
gamma = (0:10:100)
iters = [NaN,5, 5, 5, 5, 5, 5,...
         5, 5, 4, 4, 4, 4];
mask_max = [ 0.0, 0.625, 1.241, 2.491, 3.732, 4.965, 6.224,...
             7.439, 8.628, 9.805,10.883,12.543,...
            13.477,14.221,15.925,16.660,18.216,...
            19.846,21.046,21.170,22.517,23.702];
mask_min = [ 0.0, 0.556, 1.111, 2.222, 3.359, 4.497, 5.608,...
             6.658, 7.804, 9.097,10.113,11.480,...
            12.459,12.589,13.665,14.762,16.573,...
            18.214,19.312,19.641,20.597,21.656];
mask_avg = [ 0.0, 0.624, 1.238, 2.462, 3.689, 4.909, 6.122,...
             7.339, 8.513, 9.692,10.717,11.887,...
            12.951,13.790,14.975,16.135,17.520,...
            18.941,19.999,20.647,21.800,22.722];
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
xlabel('Max. % Coverage');
legend('max. coverage','min. coverage','avg. coverage','convergence')
ylim([0, 100])
hold on
%plot(fliplr(gamma),max_per,'--*',fliplr(gamma),min_per,'--x',fliplr(gamma),avg_per,'-o')
set(gca,'fontsize', fontsize);
set(gca, 'XDir','reverse')