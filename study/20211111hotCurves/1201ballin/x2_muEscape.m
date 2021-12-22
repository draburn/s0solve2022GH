clear;
thisFile = "x2_muEscape";
commondefs;
numFigs = 0;
startTime = time();

x = linspace( 0.8, 1.6, 1001 );
f = 10 - x - 0.1*x.^4;
g = 0.5*(x.^2);
%
%mu = [ 1.38, 1.39, 1.34, 1.36, 1.38, 1.40 ]
mu = linspace( 1.392, 1.393, 6 );
for m=1:max(size(mu))
	h = f+mu(m)*g;
	n = 1;
	while (n<max(size(x)))
		if ( h(n+1) > h(n) )
			break;
		end
		n++;
	end
	nCrit(m) = n;
	hCrit(m) = h(n);
end

%
plot( ...
  x, f+mu(1)*g, '-', x(nCrit(1)), hCrit(1), 'x', 'markersize', 20, 'linewidth', 3, ...
  x, f+mu(2)*g, '-', x(nCrit(2)), hCrit(2), 'x', 'markersize', 20, 'linewidth', 3, ...
  x, f+mu(3)*g, '-', x(nCrit(3)), hCrit(3), 'x', 'markersize', 20, 'linewidth', 3, ...
  x, f+mu(4)*g, '-', x(nCrit(4)), hCrit(4), 'x', 'markersize', 20, 'linewidth', 3, ...
  x, f+mu(5)*g, '-', x(nCrit(5)), hCrit(5), 'x', 'markersize', 20, 'linewidth', 3, ...
  x, f+mu(6)*g, '-', x(nCrit(6)), hCrit(6), 'x', 'markersize', 20, 'linewidth', 3 );
grid on;
