% Tarefa 1 - Alinhar histograma com curva
figure();
hold on;
printf("\n==== Tarefa 1 ===\n");
k = ceil(1 + 3.322 * log(n))
[nn_hist xx_hist] = hist(xx, k, 1);
hist(xx, k, 1)


% Tarefa 2
xx_fit = -7:0.01:3;
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate));
yy_hat = normpdf(xx_hist, mu_estimate, sigma_estimate);
RSS = sumsq(yy_hat - nn_hist)
