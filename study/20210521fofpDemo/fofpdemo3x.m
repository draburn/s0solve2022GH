clear;
% Say we have pts a, b, c, then d, and *then* the ptwise min in g.
% The local min may or may not be past d, but it must be past c.
funch_lambda = @(gamma,beta)( ...
   log(abs(gamma)./abs(gamma-beta)) ...
 ./log(abs(gamma-beta)./abs(gamma-1.0)) );
funch_lambdaLo = @(gamma,beta)( log(1.0./(1.0-beta))./log((1.0-beta)./(gamma-1.0)) );
funch_lambdaHi = @(gamma,beta)( (1.0-1.0./gamma).*beta./(1.0-beta) );
funch_lambdaMax = @(gamma,beta)( ones(size(gamma)).*beta./(1.0-beta) );

numPts = 1000;
beta = 0.01;
%gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = linspace(beta,(1.0+beta)/2.0,numPts+2);
gam1 = gamma(2:end-1);
bet1 = beta;
beta = 0.2;
%gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = linspace(beta,(1.0+beta)/2.0,numPts+2);
gam2 = gamma(2:end-1);
bet2 = beta;
beta = 0.5;
%gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = linspace(beta,(1.0+beta)/2.0,numPts+2);
gam3 = gamma(2:end-1);
bet3 = beta;
beta = 0.8;
%gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = linspace(beta,(1.0+beta)/2.0,numPts+2);
gam4 = gamma(2:end-1);
bet4 = beta;
beta = 0.95;
%gamma = linspace(beta/2.0,(1.0+beta)/2.0,numPts+2);
gamma = linspace(beta,(1.0+beta)/2.0,numPts+2);
gam5 = gamma(2:end-1);
bet5 = beta;

plot( ...
  gam1, funch_lambda(gam1,bet1), 'o-', ...
  gam2, funch_lambda(gam2,bet2), 'o-', ...
  gam3, funch_lambda(gam3,bet3), 'o-', ...
  gam4, funch_lambda(gam4,bet4), 'o-', ...
  gam5, funch_lambda(gam5,bet5), 'o-' );
grid on;
axis([0.0, 1.0, -10.0, 0.0 ]);
return;
%
beta = 0.9;
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
%f1 += (f1>10.0).*( 10.0 - f1 ) + (f1<-10.0).*( -10.0 - f1 );
plot( ...
  gamma, f1, 'o-' );
grid on;
%axis([0,1,-100,10]);
return;
gamma = linspace( 1.0, 20.0, 1000 );
beta = 0.01;
plot( ...
  gamma, funch_lambda(gamma,beta),    'o-', 'color', [0.0 0.0 0.6], ...
  gamma, funch_lambdaLo(gamma,beta),  '-',  'color', [0.6 0.6 0.9], ...
  gamma, funch_lambdaHi(gamma,beta),  '-',  'color', [0.0 0.0 1.0], ...
  gamma, funch_lambdaMax(gamma,beta), '-',  'color', [0.0 0.0 0.4] );
grid on;
axis([min(gamma), max(gamma), 0.0, 1.2*beta./(1.0-beta)] );