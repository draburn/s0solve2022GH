clear;
commondefs;
thisFile = "studyLoggyDeriv";
numFigs = 0;

caseNum = 3;
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
	bigF0 = 1e-6*abs(randn());
	bigF1 = randn();
	bigP = 2.0+abs(randn());
	funchF = @(x)( bigF0 + bigF1*( abs(x-bigX).^bigP ) );
case 2
	bigX = 0.0;
	funchF = @(x)( exp(-1./(x.^2)) );
case 3
	bigX = 0.0;
	bigF0 = 0.001;
	bigF1L = 1.0;
	bigF1R = 3.0;
	bigPL = 5.5;
	bigPR = 0.6;
	funchFL = @(x)( bigF0 + bigF1L*( abs(x-bigX).^bigPL ) );
	funchFR = @(x)( bigF0 + bigF1R*( abs(x-bigX).^bigPR ) );
	funchF = @(x)( funchFL(x) + (x>bigX).*( funchFR(x) - funchFL(x) ) );
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
numPts = 1000;
%x = bigX+0.3*linspace(-1,1,numPts)';
%x = bigX+0.5*linspace(-1,1,numPts)';
x = bigX+1.5*linspace(-1,1,numPts)';
f = funchF(x);

numFigs++; figure(numFigs);
plot( x, f, 'o-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );

for n=1:5
	dx = 10.0^(n-5);
	dxVals(n) = dx;
	xl = x-dx;
	xr = x+dx;
	fl = funchF(xl);
	fr = funchF(xr);
	df = (fr-fl)./(xr-xl);
	ddf = 2.0*(fr+fl-2.0*f)./( (xr-xl).^2 );
	neo_df(:,n) = df;
	neo_ddf(:,n) = ddf;
	fodf(:,n) = f./df;
	dfoddf(:,n) = df./ddf;
	dfoabsddf(:,n) = df./abs(ddf);
	%foddf(:,n) = f./ddf;
	%
	twoPtStencil_f = (fl+fr)/2.0;
	twoPtStencil_fodf(:,n) = twoPtStencil_f./df;
end

numFigs++; figure(numFigs);
plot( ...
  x, neo_df(:,1), 'o-', ...
  x, neo_df(:,2), 'o-', ...
  x, neo_df(:,3), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f'" );
title( "f' vs x" );

numFigs++; figure(numFigs);
plot( ...
  x, neo_ddf(:,1), 'o-', ...
  x, neo_ddf(:,2), 'o-', ...
  x, neo_ddf(:,3), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f''" );
title( "f'' vs x" );

numFigs++; figure(numFigs);
plot( ...
  x, fodf(:,1), 'o-', ...
  x, fodf(:,2), 'o-', ...
  x, fodf(:,3), 'o-', ...
  x, fodf(:,4), 'o-', ...
  x, fodf(:,5), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f/f'" );
title( "f/f' vs x" );

dfof = sign(fodf)./( eps^0.5 + abs(fodf) );
numFigs++; figure(numFigs);
plot( ...
  x, dfof(:,1), 'o-', ...
  x, dfof(:,2), 'o-', ...
  x, dfof(:,3), 'o-', ...
  x, dfof(:,4), 'o-', ...
  x, dfof(:,5), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f'/f" );
title( "f'/f vs x" );

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
xlabel( "x" );
ylabel( "f/f' mod" );
title( "f/f' mod vs x" );

numFigs++; figure(numFigs);
plot( ...
  x, dfoddf(:,1), 'o-', ...
  x, dfoddf(:,2), 'o-', ...
  x, dfoddf(:,3), 'o-', ...
  x, dfoddf(:,4), 'o-', ...
  x, dfoddf(:,5), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f'/f''" );
title( "f'/f'' vs x" );

numFigs++; figure(numFigs);
plot( ...
  x, dfoabsddf(:,1), 'o-', ...
  x, dfoabsddf(:,2), 'o-', ...
  x, dfoabsddf(:,3), 'o-', ...
  x, dfoabsddf(:,4), 'o-', ...
  x, dfoabsddf(:,5), 'o-' );
grid on;
xlabel( "x" );
ylabel( "f'/|f''|" );
title( "f'/|f''| vs x" );

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