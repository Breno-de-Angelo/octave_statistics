% Seção 1
printf("==================\n");
printf("     Seção 1\n");
printf("==================\n");


% Tarefa 1
pkg load statistics;


% Tarefa 2
printf("\n\n\n==== Tarefa 2 ===\n");
FDP_Normal = @(x, mu, sigma) 1/sqrt(2.* pi .* sigma.^2).* exp(-(x-mu).^2/(2.*sigma.^2))


% Tarefa 3
printf("\n\n\n==== Tarefa 3 ===\n");
printf("FDP_Normal(0, 0, 1) - normpdf(0, 0, 1) = %f\n", FDP_Normal(0, 0, 1) - normpdf(0, 0, 1));
printf("FDP_Normal(0, 1, 2) - normpdf(0, 1, 2) = %f\n", FDP_Normal(0, 1, 2) - normpdf(0, 1, 2));
printf("Gerando a, b, c aleatórios entre 0 e 1 ...\n");
a = rand();
b = rand();
c = rand();
printf("a = %f\nb = %f\nc = %f\n", a, b, c);
printf("FDP_Normal(a, b, c) - normpdf(a, b, c) = %f\n", FDP_Normal(a, b, c) - normpdf(a, b, c));


% Tarefa 4
xx = -5:0.01:5;
hold on;
plot(xx, normpdf(xx, 0, 1));
plot(xx, normpdf(xx, 0, sqrt(0.2)));
plot(xx, normpdf(xx, 0, sqrt(5.0)));
plot(xx, normpdf(xx, -2, sqrt(0.5)));


% Tarefa 5
printf("\n\n\n==== Tarefa 5 ===\n");
FDA_Normal = @(x, mu, sigma) 1/2 * (1 + erf((x-mu)/(sigma * sqrt(2))))
printf("FDA_Normal(0, 0, 1) - normcdf(0, 0, 1) = %f\n", FDA_Normal(0, 0, 1) - normcdf(0, 0, 1));
printf("FDA_Normal(0, 1, 2) - normcdf(0, 1, 2) = %f\n", FDA_Normal(0, 1, 2) - normcdf(0, 1, 2));
printf("Gerando a, b, c aleatórios entre 0 e 1 ...\n");
a = rand();
b = rand();
c = rand();
printf("a = %f\nb = %f\nc = %f\n", a, b, c);
printf("FDA_Normal(a, b, c) - normcdf(a, b, c) = %f\n", FDA_Normal(a, b, c) - normcdf(a, b, c));

figure()
hold on
plot(xx, normcdf(xx, 0, 1));
plot(xx, normcdf(xx, 0, sqrt(0.2)));
plot(xx, normcdf(xx, 0, sqrt(5.0)));
plot(xx, normcdf(xx, -2, sqrt(0.5)));


% Tarefa 6
printf("\n\n\n==== Tarefa 6 ===\n");
printf("FDA = 0.99  ==>   1/2 * (1 + erf((x-mu)/(sigma * sqrt(2)))) = 0.99   ==>\n");
printf("erf((x-mu)/(sigma * sqrt(2))) = (0.99 * 2) - 1   ==>   x = erfinv((0.99 * 2) - 1) * sigma * sqrt(2) + mu\n");
printf("x = erfinv(0.98) * sqrt(0.5) * sqrt(2) + (-2)\n");
printf("x = erfinv(49/50) - 2\n");
x = erfinv(49/50) - 2


% Tarefa 7



% Tarefa 8
semente = 2021;
randn('seed', semente);
n = 1000;
mu = -2;
sigma = 0.7;
xx = mu + sigma .* randn(1, n);


% Tarefa 9
printf("\n\n\n==== Tarefa 9 ===\n");
mu_estimate = 1/n * sum(xx)
sigma_estimate = sqrt(1/(n-1)*sumsq(xx - mu_estimate))
printf("|mu_estimate - mu| = %f\n",  abs(mu_estimate - mu));
printf("|sigma_estimate - sigma| = %f\n",  abs(sigma_estimate - sigma));


% Tarefa 10




% Seção 2
printf("\n\n\n");
printf("==================\n");
printf("     Seção 2\n");
printf("==================\n");


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




% Seção 3
printf("\n\n\n");
printf("==================\n");
printf("     Seção 3\n");
printf("==================\n");


% Tarefa 1
printf("\n==== Tarefa 1 ===\n");
FDP_Weibull = @(x, lambda, k) k/lambda .* (x./lambda)^(k-1) .* exp(-(x./lambda).^k)


% Tarefa 2
printf("\n\n\n==== Tarefa 2 ===\n");
printf("FDP_Weibull(1, 0.5, 1) - wblpdf(1, 0.5, 1) = %f\n", FDP_Weibull(1, 0.5, 1) - wblpdf(1, 0.5, 1));
printf("FDP_Weibull(1, 1.5, 2) - wblpdf(1, 1.5, 2) = %f\n", FDP_Weibull(1, 1.5, 2) - wblpdf(1, 1.5, 2));
printf("Gerando a, b, c aleatórios entre 0 e 1 ...\n");
a = rand();
b = rand();
c = rand();
printf("a = %f\nb = %f\nc = %f\n", a, b, c);
printf("FDP_Weibull(a, b, c) - wblpdf(a, b, c) = %f\n", FDP_Weibull(a, b, c) - wblpdf(a, b, c));


% Tarefa 3
figure();
hold on;
xx = 0:0.01:2.5;
plot(xx, wblpdf(xx, 1, 0.5));
plot(xx, wblpdf(xx, 1, 1));
plot(xx, wblpdf(xx, 1, 1.5));
plot(xx, wblpdf(xx, 1, 5));


% Tarefa 4
printf("\n\n\n==== Tarefa 4 ===\n");
FDA_Weibull = @(x, lambda, k) 1 - exp(-(x ./ lambda).^k)

printf("FDA_Weibull(1, 0.5, 1) - wblcdf(1, 0.5, 1) = %f\n", FDA_Weibull(1, 0.5, 1) - wblcdf(1, 0.5, 1));
printf("FDA_Weibull(1, 1.5, 2) - wblcdf(1, 1.5, 2) = %f\n", FDA_Weibull(1, 1.5, 2) - wblcdf(1, 1.5, 2));
printf("Gerando a, b, c aleatórios entre 0 e 1 ...\n");
a = rand();
b = rand();
c = rand();
printf("a = %f\nb = %f\nc = %f\n", a, b, c);
printf("FDA_Weibull(a, b, c) - wblcdf(a, b, c) = %f\n", FDA_Weibull(a, b, c) - wblcdf(a, b, c));

figure();
hold on;
xx = 0:0.01:2.5;
plot(xx, wblcdf(xx, 1, 0.5));
plot(xx, wblcdf(xx, 1, 1));
plot(xx, wblcdf(xx, 1, 1.5));
plot(xx, wblcdf(xx, 1, 5));


% Tarefa 5
xx = wblrnd(1, 1, [1 1000]);


% Tarefa 6 - Usar wblfit
printf("\n\n\n==== Tarefa 6 ===\n");
phi = @(xx, k) (mean(xx .^ k .* log(xx))/mean(xx .^ k) - mean(log(xx))) ^ (-1);
k_estimate = 1;
erro = inf;
while (erro > 10^(-14))
	k_old = k_estimate;
	k_estimate = phi(xx, k_estimate);
	erro = abs((k_estimate - k_old)/k_old);
endwhile
k
lambda_estimate = mean(xx .^ k_estimate) ^ (1/k_estimate)


% Tarefa 7 - Corrigir
figure();
hold on;
printf("\n\n\n==== Tarefa 7 ===\n");
k = ceil(1 + 3.322 * log(n))
[nn_hist xx_hist] = hist(xx, k, 1);
hist(xx, k, 1);

xx_fit = 0:0.01:8;
plot(xx_fit, wblpdf(xx_fit, lambda_estimate, k_estimate));
% yy_hat = wblpdf(xx_hist, lambda_estimate, k_estimate);
% RSS = sumsq(yy_hat - nn_hist)


% Tarefa 8
printf("\n\n\n==== Tarefa 8 ===\n");
mu_estimate = 1/n * sum(xx)
sigma_estimate = sqrt(1/(n-1)*sumsq(xx - mu_estimate))
printf("|mu_estimate - mu| = %f\n",  abs(mu_estimate - mu));
printf("|sigma_estimate - sigma| = %f\n",  abs(sigma_estimate - sigma));
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate));
% yy_hat = normpdf(xx_fit, mu_estimate, sigma_estimate);
% RSS = sumsq(yy_hat - nn_hist)

