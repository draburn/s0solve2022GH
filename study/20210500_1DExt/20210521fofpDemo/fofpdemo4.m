clear;
numFigs = 0;
funch_lambda = @(gamma,beta)( ...
   log(abs(gamma)./abs(gamma-beta)) ...
 ./log(abs(gamma-beta)./abs(gamma-1.0)) );
funch_lambdaLo = @(gamma,beta)( log(1.0./(1.0-beta))./log((1.0-beta)./(gamma-1.0)) );
funch_lambdaHi = @(gamma,beta)( (1.0-1.0./gamma).*beta./(1.0-beta) );
funch_lambdaMax = @(gamma,beta)( ones(size(gamma)).*beta./(1.0-beta) );
%

if (1)
	xa = 0.0
	xc = 1.0
	xb = xa + 0.05*(xc-xa);
	ga = e
	gb = 1.0
	gc = e^5
else
	setprngstates();
	xa = randn
	xc = xa+abs(randn)
	gb = abs(randn)
	%
	xb = xa + 0.05*(xc-xa);
	ga = gb*e
	gc = gb*e^5
end
bigLambda = log(ga/gb)/log(gb/gc)
beta = (xb-xa)/(xc-xa)

numPts = 1000;
gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = gamma(2:end-1);
f1 = funch_lambda(gamma,beta);
f1Sort = sort(f1);
f1PreHi = f1Sort(round(0.9*numPts));
f1PreLo = f1Sort(round(0.1*numPts));
f1Lo = f1PreLo - (f1PreHi-f1PreLo);
f1Hi = f1PreHi + (f1PreHi-f1PreLo);
f1 = cap(f1,f1Lo,f1Hi);

gamX =[0.03446, 0.1342, 0.3509]
numFigs++; figure(numFigs);
plot( ...
  gamma, f1, 'o-', ...
  gamma, bigLambda+0*gamma, 'k-', ...
  gamX(1), bigLambda, 'x', 'linewidth', 3, 'markersize', 20, ...
  gamX(2), bigLambda, '+', 'linewidth', 3, 'markersize', 22, ...
  gamX(3), bigLambda, '*', 'linewidth', 3, 'markersize', 25 );
grid on;

p = log(ga/gc)./log(gamX./(1.0-gamX))
bigX = xa + gamX*(xc-xa);
bigG = ga ./ (( gamX*(xc-xa) ).^p)

vecBigX = [ xa; xb; xc ];
matBigX = [ ones(3,1), vecBigX, vecBigX.^2 ];
vecG = [ ga; gb; gc ];
vecC = matBigX \ vecG;
c0 = vecC(1);
c1 = vecC(2);
c2 = vecC(3);

x = linspace(xa,xc,1000);
numFigs++; figure(numFigs);
plot( ...
  x, bigG(1)*abs(x-bigX(1)).^p(1), '-', 'linewidth', 2, ...
  x, bigG(2)*abs(x-bigX(2)).^p(2), '-', 'linewidth', 2, ...
  x, bigG(3)*abs(x-bigX(3)).^p(3), '-', 'linewidth', 2, ...
  x, c0 + (c1*x) + (c2*(x.^2)), '-', 'linewidth', 2, ...
  bigX(1), 0.0, 'x', 'markersize', 15, 'linewidth', 2, ...
  bigX(2), 0.0, 'x', 'markersize', 15, 'linewidth', 2, ...
  bigX(3), 0.0, 'x', 'markersize', 15, 'linewidth', 2, ...
  xa, ga, 'k*', 'markersize', 20, 'linewidth', 3, ...
  xb, gb, 'k*', 'markersize', 20, 'linewidth', 3, ...
  xc, gc, 'k*', 'markersize', 20, 'linewidth', 3 );
grid on;
