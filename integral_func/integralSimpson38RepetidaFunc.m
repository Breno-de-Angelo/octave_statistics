%%
%% Integração numérica usando a regra 3/8 de Simpson repetida
%%
%% Input: a função a ser integrada, os limites inferior 'a' e superior 'b' de integração
%%        o número 'n' de subdivisões (=número de nós menos 1).
%%
%% Output: Integral numérica
%% Restrições: número 'n' de subdivisões tem que ser divisível por três
%%
function ISR38 = integralSimpson38RepetidaFunc( func, a, b, n, verbose  )

	%% Decide se número é múltiplo de três: Uso de uma função anônima
	% http://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html
	ismult3 = @(x) (x - 3*floor(x/3) == 0);

	if ~ismult3( n )
%		fprintf('Numero de subdivisoes n=%d. \n', n); % wait();
		ISR38 = NaN;
		return;
	end
	n3 = n/3;
	h = (b-a)/n;
	nos = a:h:b;
	if verbose
		fprintf('Integracao pela Regra 3/8 de Simpson ');
		if n > 4 fprintf('repetida '), end;
		fprintf('no intervalo [%.2f,%.2f] com %d nos\n', a, b, n+1 );
		fprintf('a=%.2f, b=%.2f, n=%d, ===> h=%.2f\n', a, b, n, h);
		printTabXY( nos, 'Nos de interpolacao', func(nos), 'f(nos)', '%9.4f', 'x, f(x)' );
		xidx = 0:1:n;
		midx = ones(n+1, 1);
	end

	limits = func(a) + func(b);
	x = a + h;
	sum2 = 0.0;
	sum3 = 0.0;
	for i=1:n-1
		if ~ismult3(i)
			sum3 = sum3 + func(x); midx(i+1) = 3;
		else
			sum2 = sum2 + func(x); midx(i+1) = 2;
		end
		x = x + h;
	end
	ISR38 = 3*h/8 * (limits + 2 * sum2 + 3 * sum3);

    if ~verbose
      return;
    end
	printTabXY( xidx, 'Indice do no', midx, 'Multiplicador', '%3d', 'Multiplicadores dos nos' );
	% Explicitação didática
	fprintf('I_SR38_%d = 3h/8[ f(x0)+f(x%d) + 3{ ', n, n );
	for i=1:n3
		fprintf('f(x%d)+f(x%d)', 3*i-2, 3*i-1 ); if i==n3 fprintf(' }'); else fprintf('+'); end;
	end
	if n3 > 1 fprintf(' + 2{ '), end;
	for i=1:n3-1
		fprintf('f(x%d)', 3*i ); if i==n3-1 fprintf(' }'); else fprintf('+'); end;
	end
	fprintf(' ]\n');

	fprintf('       = 3*%.2f/8[ f(%.2f)+f(%.2f) + 3{ ', h, a, b );

	for i=1:n3
		fprintf('f(%.2f)+f(%.2f)', a+h*(3*i-2), a+h*(3*i-1) ); if i==n3 fprintf(' }'); else fprintf('+'); end;
	end
	if n3 > 1 fprintf(' + 2{ '), end;
	for i=1:n3-1
		fprintf('f(%.2f)', a+h*3*i ); if i==n3-1 fprintf(' }'); else fprintf('+'); end;
	end
	fprintf(' ]\n');

	fprintf('       = %.3f[ %.3f+%.3f + 3{ ', 3*h/8, func(a), func(b) );
	for i=1:n3
		fprintf('%.2f+%.2f', func(a+h*(3*i-2)), func(a+h*(3*i-1)) ); if i==n3 fprintf(' }'); else fprintf('+'); end;
	end
	if n3 > 1 fprintf(' + 2{ '), end;
	for i=1:n3-1
		fprintf('%.3f', func(a+h*3*i) ); if i==n3-1 fprintf(' }'); else fprintf('+'); end;
	end
	fprintf(' ] = %.10f\n', ISR38);

	numcasas = 6;
	fprintf('       = %s [ %s + %s + 3{ ', dec2fracstr(3*h/8,numcasas),...
			dec2fracstr(func(a),numcasas),...
			dec2fracstr(func(b),numcasas) );

	for i=1:n3
		fprintf('%s+%s', dec2fracstr( func(a+h*(3*i-2)), numcasas), dec2fracstr( func(a+h*(3*i-1)), numcasas) ); if i==n3 fprintf(' }'); else fprintf('+'); end;
	end
	if n3 > 1 fprintf(' + 2{ '), end;
	for i=1:n3-1
		fprintf('%s', dec2fracstr(func(a+h*3*i), numcasas) ); if i==n3-1 fprintf(' }'); else fprintf('+'); end;
	end
	fprintf(' ] = '); printdecandfrac( ISR38, true );
end

