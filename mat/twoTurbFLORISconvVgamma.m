close()
fontsize = 14;
gamma = [0.0,0.01,0.015,0.025,0.05,0.1:0.1:1.0];
max_per = [100,98.75,94.31,41.53,36.25,31.50, 26.25, 20.97, 18.10,14.02,8.89,7.08,3.47,2.22,0];
min_per = [100,38.19,34.86,33.06,29.44,25.83,20.69,14.44,10.97,6.25,2.08,0.97,0.28,0.14,0];
avg_per = [100,62.25,48.72,36.79,33.16,28.80,23.52,18.19,14.13,9.82,6.52,3.70,1.74,0.32,0];
converg = [5,5,5,7,8,9,11,14,17,23,34,61,135,400,100000000];
yyaxis left
plot(gamma,max_per,'--*',gamma,min_per,'--x',gamma,avg_per,'-o')
ylabel('% coverage')
xlabel('\gamma')
yyaxis right
plot(gamma,converg,'-.o')
ylim([0,400])
legend('max. coverage','min. coverage','avg. coverage','convergence')
ylabel("steps to convergence")
set(gca,'fontsize', fontsize);