clear;
commondefs;
thisFile = "studyLoggyDeriv";
numFigs = 0;

caseNum = 0;
switch (caseNum)
case 0
	bigX = 0.0;
	bigF0 = 0.0;
	bigF1 = 1.0;
	bigP = 10.0;
	funchF = @(x)( bigF0 + bigF1*( abs(x-bigX).^bigP ) );
case 1
	setprngstates(0);
	bigX = randn()*exp(abs(randn()));
	bigF0 = randn();
	bigF1 = randn();
	bigP = abs(randn());
	funchF = @(x)( bigF0 + bigF1*( abs(x-bigX).^bigP ) );
case 2
	funchF = @(x)( exp(-1./(x.^2)) );
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
numPts = 1000;
x = 2*linspace(-1,1,numPts)';
f = funchF(x);

numFigs++; figure(numFigs);
plot( x, f, 'o-' );
grid on;

for n=1:5
	dx = 10.0^(n-5);
	xl = x-dx;
	xr = x+dx;
	fl = funchF(xl);
	fr = funchF(xr);
	df = (fr-fl)./(xr-xl);
	ddf = 2.0*(fr+fl-2.0*f)./( (xr-xl).^2 );
	fodf(:,n) = f./df;
	dfoddf(:,n) = df./ddf;
	%foddf(:,n) = f./ddf;
	%
	twoPtStencil_f = (fl+fr)/2.0;
	twoPtStencil_fodf(:,n) = twoPtStencil_f./df;
end

numFigs++; figure(numFigs);
plot( ...
  x, fodf(:,1), 'o-', ...
  x, fodf(:,2), 'o-', ...
  x, fodf(:,3), 'o-', ...
  x, fodf(:,4), 'o-', ...
  x, fodf(:,5), 'o-' );
grid on;

for n=1:5
	fodf_mod(:,n) = fodf(:,n).*( 1.0 + 0.001*abs(x).^(-2) );
end
numFigs++; figure(numFigs);
plot( ...
  x, fodf_mod(:,1), 'o-', ...
  x, fodf_mod(:,2), 'o-', ...
  x, fodf_mod(:,3), 'o-', ...
  x, fodf_mod(:,4), 'o-', ...
  x, fodf_mod(:,5), 'o-' );
grid on;

numFigs++; figure(numFigs);
plot( ...
  x, dfoddf(:,1), 'o-', ...
  x, dfoddf(:,2), 'o-', ...
  x, dfoddf(:,3), 'o-', ...
  x, dfoddf(:,4), 'o-', ...
  x, dfoddf(:,5), 'o-' );
grid on;

%numFigs++; figure(numFigs);
%plot( ...
%  x, twoPtStencil_fodf(:,1), 'o-', ...
%  x, twoPtStencil_fodf(:,2), 'o-', ...
%  x, twoPtStencil_fodf(:,3), 'o-', ...
%  x, twoPtStencil_fodf(:,4), 'o-', ...
%  x, twoPtStencil_fodf(:,5), 'o-' );
%grid on;

%numFigs++; figure(numFigs);
%plot( ...
%  x, foddf(:,1), 'o-', ...
%  x, foddf(:,2), 'o-', ...
%  x, foddf(:,3), 'o-', ...
%  x, foddf(:,4), 'o-', ...
%  x, foddf(:,5), 'o-' );
%grid on;

return;

dx = 1.0;
xm2 = x-(2*dx);
xm1 = x-dx;
xp1 = x+dx;
xp2 = x+(2*dx);
fm2 = funchF(xm2);
fm1 = funchF(xm1);
fp1 = funchF(xp1);
fp2 = funchF(xp2);
quad_df = zeros(size(f));
cubi_df = zeros(size(f));
for n=1:numPts
	vecF = [ fm2(n); fm1(n); f(n); fp1(n); fp2(n) ];
	vecX = [ xm2(n); xm1(n); x(n); xp1(n); xp2(n) ];
	matX = [ ones(size(vecX)), vecX, vecX.^2 ];
	vecC = matX \ vecF;
	quad_df(n) = vecC(2) + 2.0*vecC(3)*x(n);
	%
	%matX = [ ones(size(vecX)), vecX, vecX.^2, vecX.^3 ];
	%vecC = matX \ vecF;
end
quad_fodf = f./quad_df;


numFigs++; figure(numFigs);
plot( ...
  x, fodf(:,5), 'o-', ...
  x, quad_fodf, 'o-' );
grid on;