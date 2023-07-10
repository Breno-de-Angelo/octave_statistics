addpath(['./' 'equation_root'], ['./' 'integral_func'])


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

% Adicionar legenda
legend('\mu = 0, \sigma^2 =1', '\mu = 0, \sigma^2 =0.2', '\mu = 0,\sigma^2 =5.0', '\mu = -2,\sigma^2 = 0.5');

% Nome nos eixos e título
xlabel('x');
ylabel(' \phi_{\mu,\sigma^2}(x)');
title('Distribuições Normais');

hold off;


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

% Adicionar legenda
legend('\mu = 0, \sigma^2 =1', '\mu = 0, \sigma^2 =0.2', '\mu = 0,\sigma^2 =5.0', '\mu = -2,\sigma^2 = 0.5', 'Location', 'southeast');

% Nome nos eixos e título
xlabel('x');
ylabel(' \Phi_{\mu,\sigma^2}(x)');
title('Função de Distribuição Acumulada');



% Tarefa 6
printf("\n\n\n==== Tarefa 6 ===\n");
printf("FDA = 0.99  ==>   1/2 * (1 + erf((x-mu)/(sigma * sqrt(2)))) = 0.99   ==>\n");
printf("erf((x-mu)/(sigma * sqrt(2))) = (0.99 * 2) - 1   ==>   x = erfinv((0.99 * 2) - 1) * sigma * sqrt(2) + mu\n");
printf("x = erfinv(0.98) * sqrt(0.5) * sqrt(2) + (-2)\n");
printf("x = erfinv(49/50) - 2\n");
x = erfinv(49/50) - 2


% Tarefa 7
mu = -2;
sigma2 = 0.5;
sigma = sqrt(sigma2);
func = @(x) FDP_Normal(x, mu, sigma);
a = mu - 10 * sigma;
b = -2 + erfinv(49/50);
erros = zeros([12 4]);
fprintf('------------------------------------------------------------------------------------------\n');
fprintf('|   n   |   Trapezios   |   Simpson   |   Simpson 3/8   |   Quadratura Gaussiana   |\n');
fprintf('------------------------------------------------------------------------------------------\n');

for n = 1:12
    fprintf('|   %2d  |', n);

	integral_correct = normcdf(b, mu, sigma) - normcdf(a, mu, sigma);

    integral_result = integralTrapeziosRepetidaFunc(func, a, b, n, 0);
    error_trapezoid = abs(integral_correct - integral_result);
    %fprintf('   %12.6f   |', integral_result);
    fprintf('   %12.6f   |', error_trapezoid);
    erros(n, 1) = error_trapezoid;

    integral_result = integralSimpsonRepetidaFunc(func, a, b, n, 0);
    error_simpson = abs(integral_correct - integral_result);
    %fprintf('   %12.6f   |', integral_result);
    fprintf('   %12.6f   |', error_simpson);
    erros(n, 2) = error_simpson;

    integral_result = integralSimpson38RepetidaFunc(func, a, b, n, 0);
    error_simpson38 = abs(integral_correct - integral_result);
    %fprintf('   %15.6f   |', integral_result);
    fprintf('   %15.6f   |', error_simpson38);
    erros(n, 3) = error_simpson38;

    C = coefGaussLegendre(n + 1);
    [T, A] = tabelaAbcissasPesosGaussLegendre(C);
    integral_result = integralGaussLegendreFunc(func, a, b, n, T, A, 0);
    error_gauss = abs(integral_correct - integral_result);
    %fprintf('   %22.6f   |', integral_result);
    fprintf('   %22.6f   |', error_gauss);
    erros(n, 4) = error_gauss;

    fprintf('\n');
end
figure()
bar(1:12, erros);
legend('Trapezios', 'Simpson 1/3', 'Simpson 3/8', 'Quadratura Gaussiana');
xlabel('Quantidade de Subdivisões');
ylabel('Erro');
title('Evolução do Erro em Relação à Quantidade de Subdivisões');
fprintf('------------------------------------------------------------------------------------------\n');


% Tarefa 8
semente = 2021;
randn('seed', semente);
n = 1000;
xx = mu + sigma .* randn(1, n);


% Tarefa 9
printf("\n\n\n==== Tarefa 9 ===\n");
mu_estimate = 1/n * sum(xx)
sigma_estimate = sqrt(1/(n-1)*sumsq(xx - mu_estimate))
printf("|mu_estimate - mu| = %f\n",  abs(mu_estimate - mu));
printf("|sigma_estimate - sigma| = %f\n",  abs(sigma_estimate - sigma));


% Tarefa 10
printf("\n\n\n==== Tarefa 9 ===\n");
func = @(x) FDA_Normal(x, mu, sigma) - 0.99;
raizBisecPosFalsa( 0, func, -1, 0, 0, 10000)
solucao_analitica = mu + sigma * sqrt(2) * erfinv(2 * 0.99 - 1)



% Seção 2
printf("\n\n\n");
printf("==================\n");
printf("     Seção 2\n");
printf("==================\n");


% Tarefa 1
printf("\n\n\n==== Tarefa 1 ===\n");
figure();
hold on;
k = ceil(1 + 3.322 * log(n))
[nn_hist xx_hist] = hist(xx, k);
delta_x = xx_hist(2) - xx_hist(1);
[yy_hist xx_hist] = hist(xx, k, 1/delta_x);
hist(xx, k, 1/delta_x)
xlabel('Distribuicao Normal');
ylabel('Densidade de Probabilidade');
title('Histograma');


% Tarefa 2
printf("\n\n\n==== Tarefa 2 ===\n");
xx_fit = -7:0.01:3;
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate));
yy_hat = normpdf(xx_hist, mu_estimate, sigma_estimate);
RSS = sumsq(yy_hat - yy_hist)



% Seção 3
printf("\n\n\n");
printf("==================\n");
printf("     Seção 3\n");
printf("==================\n");


% Tarefa 1
printf("\n\n\n==== Tarefa 1 ===\n");
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

% Adicionar legenda
legend('\lambda=1,k=0,5', '\lambda=1,k=1','\lambda=1,k=1,5', '\lambda=1,k=5');

% Nome nos eixos e título
xlabel('s');
ylabel('P');
title('Distribuicao de Wibull');



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

% Adicionar legenda
legend('\lambda=1,k=0,5', '\lambda=1,k=1','\lambda=1,k=1,5', '\lambda=1,k=5', 'Location', 'southeast');

% Nome nos eixos e título
xlabel('s');
ylabel('P');
title('Distribuicao de Wibull FDA');


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


% Tarefa 7
figure();
hold on;
printf("\n\n\n==== Tarefa 7 ===\n");
k = ceil(1 + 3.322 * log(n));
[nn_hist xx_hist] = hist(xx, k);
delta_x = xx_hist(2) - xx_hist(1);
[yy_hist xx_hist] = hist(xx, k, 1/delta_x);
hist(xx, k, 1/delta_x);
xx_fit = 0:0.01:8;
plot(xx_fit, wblpdf(xx_fit, lambda_estimate, k_estimate));
yy_hat = wblpdf(xx_hist, lambda_estimate, k_estimate);
RSS = sumsq(yy_hat - yy_hist);





% Tarefa 8
printf("\n\n\n==== Tarefa 8 ===\n");
mu_estimate = 1/n * sum(xx)
sigma_estimate = sqrt(1/(n-1)*sumsq(xx - mu_estimate))
printf("|mu_estimate - mu| = %f\n",  abs(mu_estimate - mu));
printf("|sigma_estimate - sigma| = %f\n",  abs(sigma_estimate - sigma));
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate));
yy_hat = normpdf(xx_hist, mu_estimate, sigma_estimate);
RSS = sumsq(yy_hat - yy_hist)

% Adicionar legenda
legend('y observado','y predito pela FDP','ajuste com a distribuicao proposta');

xlabel('Intervalos');
ylabel('Frequência');
title('Estimativa da FDP com Mínimos Quadrados e Histograma');
grid on;


% Seção 4
printf("\n\n\n");
printf("==================\n");
printf("     Seção 4\n");
printf("==================\n");


% Tarefa 1
printf("\n\n\n==== Tarefa 1 ===\n");
FDP_Pareto = @(x, mu, sigma, ksi) 1/sigma * (1 + ksi * (x - mu) / sigma) .^ (-1/ksi - 1)


% Tarefa 2
xx = 0:0.01:5;
figure();
hold on;
plot(xx, FDP_Pareto(xx, 0, 1, 1),'r');
plot(xx, FDP_Pareto(xx, 0, 1, 5),'g');
plot(xx, FDP_Pareto(xx, 0, 1, 20),'b');
plot(xx, FDP_Pareto(xx, 0, 2, 1),'r--');
plot(xx, FDP_Pareto(xx, 0, 2, 5),'g--');
plot(xx, FDP_Pareto(xx, 0, 2, 20),'b--');


% Adicionar legenda
legend('\sigma = 1, \xi = 1','\sigma = 1,\xi = 5', '\sigma = 1, \xi = 20','\sigma = 2, \xi = 1','\sigma =1, \xi = 1', '\sigma = 2, \xi = 5', '\sigma = 2, \xi = 20');

xlabel('s');
ylabel('P');
title('Distribuicao Generalizada de Pareto para FDP');

% Tarefa 3
function retval = FDA_Pareto(x, mu, sigma, ksi)
	if (ksi == 0)
		retval = 1 - exp(- (x - mu) / sigma);
	else
		retval = 1 - (1 + ksi * (x - mu) / sigma) .^ (-1/ksi);
	endif
endfunction
figure();
hold on;
plot(xx, FDA_Pareto(xx, 0, 1, 1),'r');
plot(xx, FDA_Pareto(xx, 0, 1, 5),'g');
plot(xx, FDA_Pareto(xx, 0, 1, 20),'b');
plot(xx, FDA_Pareto(xx, 0, 2, 1),'r--');
plot(xx, FDA_Pareto(xx, 0, 2, 5),'g--');
plot(xx, FDA_Pareto(xx, 0, 2, 20),'b--');


% Adicionar legenda
legend('\sigma = 1, \xi = 1','\sigma = 1,\xi = 5', '\sigma = 1, \xi = 20','\sigma = 2, \xi = 1','\sigma =1, \xi = 1', '\sigma = 2, \xi = 5', '\sigma = 2, \xi = 20');

xlabel('s');
ylabel('P');
title('Distribuicao Generalizada de Pareto para FDA');

% Tarefa 4
n = 1000;
xx = gprnd(1, 1, 0, [1 n]);


% Tarefa 5
printf("\n\n\n==== Tarefa 5 ===\n");
figure();
k = ceil(1 + 3.322 * log(n))
[nn_hist xx_hist] = hist(xx, k);
delta_x = xx_hist(2) - xx_hist(1);
[yy_hist xx_hist] = hist(xx, k, 1/delta_x);
mu_estimate = 0
expected_value = 1 / n * sum(xx)
variance = 1 / (n - 1) * sum ((xx - expected_value) .^ 2)
ksi_estimate = 1/2 * (1 - (expected_value - mu_estimate) ^ 2 / variance)
sigma_estimate =  (expected_value - 0) * (1 - ksi_estimate)
xx_fit = 0:0.01:12;
yy_hat = FDP_Pareto(xx_hist, mu_estimate, sigma_estimate, ksi_estimate);
plot(xx_fit, FDP_Pareto(xx_fit, mu_estimate, sigma_estimate, ksi_estimate));
RSS = sumsq(yy_hat - yy_hist)
xlabel('s');
ylabel('P');
title('FDP sujeito a Distribuicao de Paretto Generalizada');

% PROBLEMA PRÁTICO
printf("\n\n\n");
printf("==================\n");
printf("     Problema Prático\n");
printf("==================\n");

wspd = load("vetorWSPD").wspd;
wvht = load("vetorWVHT").wvht;
wspd(wspd==0) = 0.1;
wvht(wvht==0) = 0.1;

wspd = wspd(wspd~=99);
wvht = wvht(wvht~=99);

t_wspd = [1:1:length(wspd)];
t_wvht = [1:1:length(wvht)];

figure()
plot(t_wspd, wspd);

figure()
plot(t_wvht, wvht, '.', 'MarkerSize', 0.2);


% WSPD
k_wspd = ceil(1 + 3.322 * log(length(wspd)));
mu_estimate = 1/length(wspd) * sum(wspd);
sigma_estimate = sqrt(1/(length(wspd)-1)*sumsq(wspd - mu_estimate));

figure()
[nn, xx] = hist(wspd, k_wspd);
hist(wspd, k_wspd);
delta = xx(2) - xx(1);
areaHist = sum(delta*nn);
hold on
xx_fit = xx(1):0.01:xx(end);
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate)*areaHist, '.');
plot(xx_fit, gppdf(xx_fit, gpfit(wspd)(1), gpfit(wspd)(2),0)*areaHist, '.');
hold off

% WVHT
k_wvht = ceil(1 + 3.322 * log(length(wvht)));
mu_estimate = 1/length(wvht) * sum(wvht);
sigma_estimate = sqrt(1/(length(wvht)-1)*sumsq(wvht - mu_estimate));

figure()
[nn, xx] = hist(wvht, k_wvht);
hist(wvht, k_wvht);
delta = xx(2) - xx(1);
areaHist = sum(delta*nn);
hold on
xx_fit = xx(1):0.01:xx(end);
plot(xx_fit, normpdf(xx_fit, mu_estimate, sigma_estimate)*areaHist, '.');
##plot(xx_fit, wblpdf(xx_fit, wblfit(wvht)(2), wblfit(wvht)(1))*areaHist, '.');
plot(xx_fit, gppdf(xx_fit, gpfit(wvht)(1), gpfit(wvht)(2),0)*areaHist, '.');
hold off
