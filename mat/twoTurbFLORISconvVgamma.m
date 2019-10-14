close()
fontsize = 14;
gamma = [0.0,0.01,0.015,0.025,0.05,0.1:0.1:1.0];
max_per = [100,40.0,37.08,35.56,32.5,28.19, 23.33, 18.33, 15.97,13.47,8.89,6.94,5.83,2.78,0];
min_per = [100,32.5,33.61,32.08,23.75,20.97,19.31,14.44,9.72,5.55,4.02,2.36,1.67,0.42,0];
avg_per = [100,36.97,35.56,34.0,30.03,25.97,21.25,17.01,13.17,9.95,6.39,4.65,4.23,1.08,0];
converg = [5,5,5,5,5,5,5,6,6,6,6,7,15,22,50];
yyaxis left
plot(gamma,max_per,'--*',gamma,min_per,'--x',gamma,avg_per,'-o')
ylabel('% coverage')
xlabel('\gamma')
yyaxis right
plot(gamma,converg,'-.o')
ylim([0,25])
legend('max. coverage','min. coverage','avg. coverage','convergence')
ylabel("steps to convergence")
set(gca,'fontsize', fontsize);