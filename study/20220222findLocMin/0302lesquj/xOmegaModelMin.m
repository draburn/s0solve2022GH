clear;
commondefs;
numFigs = 0;
tic();
%
%function [ vecF, matJ ] = funcFJ_easy( vecX )
%	sizeX = size(vecX,1);
%	matJ = diag((1:sizeX));
%	vecF = matJ*( vecX - (1:sizeX)' );
%endfunction
function [ vecF, matJ ] = funcFJ_nonlin( vecX )
	sizeX = size(vecX,1);
	vecXE = (1:sizeX)';
	matA = diag((1:sizeX));
	matA(2,1) = (sizeX+1);
	vecF = matA*((vecX-vecXE)).^2;
	matJ = zeros(sizeX,sizeX);
	for n=1:sizeX
	for m=1:sizeX
		matJ(n,m) = 2.0*matA(n,m)*(vecX(m)-vecXE(m));
	end
	end
endfunction
%
%
funchFJ = @(dummyX) funcFJ_nonlin(dummyX);
%%%funchFJ = @(dummyX) funcFJ_easy(dummyX);
vecX0 = zeros(2,1);
[ vecF0, matJ0 ] = funchFJ( vecX0 );
%
prm_sansCDL = [];
prm_sansCDL.useCDL = false;
[ vecXF_sansCDL, datOut_sansCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
%
prm_withCDL = [];
prm_withCDL.useCDL = true;
[ vecXF_withCDL, datOut_withCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_sansCDL.fevalCountVals, datOut_sansCDL.omegaVals, 'o-', ...
  datOut_withCDL.fevalCountVals, datOut_withCDL.omegaVals, 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "omega" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "omega vs feval count" );
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_sansCDL.fevalCountVals(2:end), sqrt(sumsq(datOut_sansCDL.vecXVals(:,2:end)-datOut_sansCDL.vecXVals(:,1:end-1),1)), 'o-', ...
  datOut_withCDL.fevalCountVals(2:end), sqrt(sumsq(datOut_withCDL.vecXVals(:,2:end)-datOut_withCDL.vecXVals(:,1:end-1),1)), 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "delta norm" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "delta norm vs feval count" );
%
numFigs++; figure(numFigs);
plot( ...
  datOut_sansCDL.fevalCountVals, datOut_sansCDL.matJVals(1,1,:), 'o-', ...
  datOut_withCDL.fevalCountVals, datOut_withCDL.matJVals(1,1,:), 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "j11" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "j11 vs feval count" );
numFigs++; figure(numFigs);
plot( ...
  datOut_sansCDL.fevalCountVals, datOut_sansCDL.matJVals(1,2,:), 'o-', ...
  datOut_withCDL.fevalCountVals, datOut_withCDL.matJVals(1,2,:), 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "j12" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "j12 vs feval count" );
numFigs++; figure(numFigs);
plot( ...
  datOut_sansCDL.fevalCountVals, datOut_sansCDL.matJVals(2,1,:), 'o-', ...
  datOut_withCDL.fevalCountVals, datOut_withCDL.matJVals(2,1,:), 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "j21" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "j21 vs feval count" );
numFigs++; figure(numFigs);
plot( ...
  datOut_sansCDL.fevalCountVals, datOut_sansCDL.matJVals(2,2,:), 'o-', ...
  datOut_withCDL.fevalCountVals, datOut_withCDL.matJVals(2,2,:), 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "j22" );
legend( "sans CDL", "with CDL", "location", "northeast" );
title( "j22 vs feval count" );
%
%
toc();
msg( __FILE__, __LINE__, "This code shows that calcDeltaLev helps make the Jacobian more accurate, which is perhaps bad here." );
msg( __FILE__, __LINE__, "In some sense, the non-{CDL+omegaModelMin} approach may just *happen* to be working well." );
msg( __FILE__, __LINE__, "This suggests that Broyden behaves poorly near a bad local min." );
msg( __FILE__, __LINE__, "Also, this may be an atypical test case." );
msg( __FILE__, __LINE__, "Let's look at a nonlinear case that has a root where the Jacobian is non-singular." );
