clear;
thisFile = "test_extFit__getLocalModel";
msg( thisFile, __LINE__, "DEPRECATED." );
bigA = 1.0;
bigB = 1.0;
bigX = 0.5;
bigP = 3.0;
funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
xVals = linspace( -2.0, 3.0, 6 );
fVals = funchF(xVals);
bigX_guess = bigX + 0.2;
bigP_guess = bigP + 0.2;
%
dat_localModel_rhoLin    = extFit__getLocalModel_rhoLin(    bigX_guess, bigP_guess, xVals, fVals );
dat_localModel_rhoSqQuad = extFit__getLocalModel_rhoSqQuad( bigX_guess, bigP_guess, xVals, fVals );
dat_localModel_omegaQuad = extFit__getLocalModel_omegaQuad( bigX_guess, bigP_guess, xVals, fVals );
%echo__rhoLin_matH = dat_localModel_rhoLin.matH
%echo__rhoSqQuad_matH = dat_localModel_rhoSqQuad.matH
%echo__omegaQuad_matH = dat_localModel_omegaQuad.matH
%
funch_omega = @(deltaVec)( 0.5*sum(extFit__getRhoVals(bigX_guess+deltaVec(1),bigP_guess+deltaVec(2),xVals,fVals).^2) );
lambdaVals = linspace(0.0,1.0,31).^2;
n = 0;
for lambda=lambdaVals
	n++;
	%
	vecDeltaVals_rhoLin_newton(:,n)     = dat_localModel_rhoLin.dat_funchDelta.newton(     lambda );
	vecDeltaVals_rhoLin_gradDir(:,n)    = dat_localModel_rhoLin.dat_funchDelta.gradDir(    lambda );
	vecDeltaVals_rhoLin_gradSclDir(:,n) = dat_localModel_rhoLin.dat_funchDelta.gradSclDir( lambda );
	vecDeltaVals_rhoLin_levenberg(:,n)  = dat_localModel_rhoLin.dat_funchDelta.levenberg(  lambda );
	vecDeltaVals_rhoLin_levMarq(:,n)    = dat_localModel_rhoLin.dat_funchDelta.levMarq(    lambda );
	%
	vecDeltaVals_rhoSqQuad_newton(:,n)     = dat_localModel_rhoSqQuad.dat_funchDelta.newton(     lambda );
	vecDeltaVals_rhoSqQuad_gradDir(:,n)    = dat_localModel_rhoSqQuad.dat_funchDelta.gradDir(    lambda );
	vecDeltaVals_rhoSqQuad_gradSclDir(:,n) = dat_localModel_rhoSqQuad.dat_funchDelta.gradSclDir( lambda );
	vecDeltaVals_rhoSqQuad_levenberg(:,n)  = dat_localModel_rhoSqQuad.dat_funchDelta.levenberg(  lambda );
	vecDeltaVals_rhoSqQuad_levMarq(:,n)    = dat_localModel_rhoSqQuad.dat_funchDelta.levMarq(    lambda );
	%
	vecDeltaVals_omegaQuad_newton(:,n)     = dat_localModel_omegaQuad.dat_funchDelta.newton(     lambda );
	vecDeltaVals_omegaQuad_gradDir(:,n)    = dat_localModel_omegaQuad.dat_funchDelta.gradDir(    lambda );
	vecDeltaVals_omegaQuad_gradSclDir(:,n) = dat_localModel_omegaQuad.dat_funchDelta.gradSclDir( lambda );
	vecDeltaVals_omegaQuad_levenberg(:,n)  = dat_localModel_omegaQuad.dat_funchDelta.levenberg(  lambda );
	vecDeltaVals_omegaQuad_levMarq(:,n)    = dat_localModel_omegaQuad.dat_funchDelta.levMarq(    lambda );
	%
	%
	omegaModelVals_rhoLin_newton(n)    = dat_localModel_rhoLin.funchOmegaModel(vecDeltaVals_rhoLin_newton(:,n));
	omegaModelVals_rhoLin_levenberg(n) = dat_localModel_rhoLin.funchOmegaModel(vecDeltaVals_rhoLin_levenberg(:,n));
	omegaModelVals_rhoLin_levMarq(n)   = dat_localModel_rhoLin.funchOmegaModel(vecDeltaVals_rhoLin_levMarq(:,n));
	%
	omegaModelVals_rhoSqQuad_newton(n)    = dat_localModel_rhoSqQuad.funchOmegaModel(vecDeltaVals_rhoSqQuad_newton(:,n));
	omegaModelVals_rhoSqQuad_levenberg(n) = dat_localModel_rhoSqQuad.funchOmegaModel(vecDeltaVals_rhoSqQuad_levenberg(:,n));
	omegaModelVals_rhoSqQuad_levMarq(n)   = dat_localModel_rhoSqQuad.funchOmegaModel(vecDeltaVals_rhoSqQuad_levMarq(:,n));
	%
	omegaModelVals_omegaQuad_newton(n)    = dat_localModel_omegaQuad.funchOmegaModel(vecDeltaVals_omegaQuad_newton(:,n));
	omegaModelVals_omegaQuad_levenberg(n) = dat_localModel_omegaQuad.funchOmegaModel(vecDeltaVals_omegaQuad_levenberg(:,n));
	omegaModelVals_omegaQuad_levMarq(n)   = dat_localModel_omegaQuad.funchOmegaModel(vecDeltaVals_omegaQuad_levMarq(:,n));
	%
	%
	omegaActualVals_rhoLin_newton(n)    = funch_omega(vecDeltaVals_rhoLin_newton(:,n));
	omegaActualVals_rhoLin_levenberg(n) = funch_omega(vecDeltaVals_rhoLin_levenberg(:,n));
	omegaActualVals_rhoLin_levMarq(n)   = funch_omega(vecDeltaVals_rhoLin_levMarq(:,n));
	%
	omegaActualVals_rhoSqQuad_newton(n)    = funch_omega(vecDeltaVals_rhoSqQuad_newton(:,n));
	omegaActualVals_rhoSqQuad_levenberg(n) = funch_omega(vecDeltaVals_rhoSqQuad_levenberg(:,n));
	omegaActualVals_rhoSqQuad_levMarq(n)   = funch_omega(vecDeltaVals_rhoSqQuad_levMarq(:,n));
	%
	omegaActualVals_omegaQuad_newton(n)    = funch_omega(vecDeltaVals_omegaQuad_newton(:,n));
	omegaActualVals_omegaQuad_levenberg(n) = funch_omega(vecDeltaVals_omegaQuad_levenberg(:,n));
	omegaActualVals_omegaQuad_levMarq(n)   = funch_omega(vecDeltaVals_omegaQuad_levMarq(:,n));
end

numFigs = 0;

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, vecDeltaVals_rhoLin_newton(1,:), '*-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(1,:), 'o-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(1,:), '^-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(1,:), '*-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(1,:), 'o-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(1,:), '^-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(1,:), '*-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(1,:), 'o-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(1,:), '^-', 'markersize', 9 );
title( "LEGEND" );
xlabel( "LEGEND" );
ylabel( "LEGEND" );
legend( ...
  "rhoLin - newt", ...
  "rhoLin - lev", ...
  "rhoLin - levMarq", ...
  "rhoSqQuad - newt", ...
  "rhoSqQuad - lev", ...
  "rhoSqQuad - levMarq", ...
  "omegaQuad - newt", ...
  "omegaQuad - lev", ...
  "omegaQuad - levMarq", ...
  "location", "northwest" );

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, vecDeltaVals_rhoLin_newton(1,:), '*-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(1,:), 'o-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(1,:), '^-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(1,:), '*-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(1,:), 'o-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(1,:), '^-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(1,:), '*-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(1,:), 'o-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(1,:), '^-', 'markersize', 9 );
grid on;
title( "deltaX vs lambda" );
xlabel( "lambda" );
ylabel( "deltaX" );

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, vecDeltaVals_rhoLin_newton(2,:), '*-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(2,:), 'o-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(2,:), '^-', 'markersize', 16, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(2,:), '*-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(2,:), 'o-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(2,:), '^-', 'markersize', 12, ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(2,:), '*-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(2,:), 'o-', 'markersize', 9, ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(2,:), '^-', 'markersize', 9 );
grid on;
title( "deltaP vs lambda" );
xlabel( "lambda" );
ylabel( "deltaP" );

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, omegaModelVals_rhoLin_newton, '*-', 'markersize', 16, ...
  lambdaVals, omegaModelVals_rhoLin_levenberg, 'o-', 'markersize', 16, ...
  lambdaVals, omegaModelVals_rhoLin_levMarq, '^-', 'markersize', 16, ...
  lambdaVals, omegaModelVals_rhoSqQuad_newton, '*-', 'markersize', 12, ...
  lambdaVals, omegaModelVals_rhoSqQuad_levenberg, 'o-', 'markersize', 12, ...
  lambdaVals, omegaModelVals_rhoSqQuad_levMarq, '^-', 'markersize', 12, ...
  lambdaVals, omegaModelVals_omegaQuad_newton, '*-', 'markersize', 10, ...
  lambdaVals, omegaModelVals_omegaQuad_levenberg, 'o-', 'markersize', 10, ...
  lambdaVals, omegaModelVals_omegaQuad_levMarq, '^-', 'markersize', 10 );
grid on;
title( "omegaModel vs lambda" );
xlabel( "lambda" );
ylabel( "omegaModel" );

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, omegaActualVals_rhoLin_newton, '*-', 'markersize', 16, ...
  lambdaVals, omegaActualVals_rhoLin_levenberg, 'o-', 'markersize', 16, ...
  lambdaVals, omegaActualVals_rhoLin_levMarq, '^-', 'markersize', 16, ...
  lambdaVals, omegaActualVals_rhoSqQuad_newton, '*-', 'markersize', 12, ...
  lambdaVals, omegaActualVals_rhoSqQuad_levenberg, 'o-', 'markersize', 12, ...
  lambdaVals, omegaActualVals_rhoSqQuad_levMarq, '^-', 'markersize', 12, ...
  lambdaVals, omegaActualVals_omegaQuad_newton, '*-', 'markersize', 10, ...
  lambdaVals, omegaActualVals_omegaQuad_levenberg, 'o-', 'markersize', 10, ...
  lambdaVals, omegaActualVals_omegaQuad_levMarq, '^-', 'markersize', 10 );
grid on;
title( "omegaActual vs lambda" );
xlabel( "lambda" );
ylabel( "omegaActual" );



numBigXVals = 21;
numBigPVals = 31;
bigXVals = linspace(bigX-1.0,bigX+1.0,numBigXVals);
bigPVals = linspace(bigP-0.5,bigP+0.5,numBigPVals);
[ mesh_bigX, mesh_bigP ] = meshgrid( bigXVals, bigPVals );
mesh_omega = zeros(numBigPVals,numBigXVals);
for ix=1:numBigXVals
for ip=1:numBigPVals
	mesh_omega(ip,ix) = 0.5*sum( extFit__getRhoVals( bigXVals(ix), bigPVals(ip), xVals, fVals ).^2 );
	mesh_omegaModel_rhoLin(ip,ix) = dat_localModel_rhoLin.funchOmegaModel([ bigXVals(ix)-bigX_guess; bigPVals(ip)-bigP_guess ]);
	mesh_omegaModel_rhoSqQuad(ip,ix) = dat_localModel_rhoSqQuad.funchOmegaModel([ bigXVals(ix)-bigX_guess; bigPVals(ip)-bigP_guess ]);
	mesh_omegaModel_omegaQuad(ip,ix) = dat_localModel_omegaQuad.funchOmegaModel([ bigXVals(ix)-bigX_guess; bigPVals(ip)-bigP_guess ]);
end
end

numFigs++; figure(numFigs);
contourf( bigXVals, bigPVals, sqrt(abs(mesh_omega)), 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "omega vs bigX, bigP" );

numFigs++; figure(numFigs);
ddx_mesh_omega = ...
 ( mesh_omega(:,2:end) - mesh_omega(:,1:end-1) ) ./ ...
 ( mesh_bigX(:,2:end) - mesh_bigX(:,1:end-1) );
contourf( cent(bigXVals), bigPVals, ddx_mesh_omega, 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "D/DX omega - rhoLin vs bigX, bigP" );

numFigs++; figure(numFigs);
ddp_mesh_omega = ...
 ( mesh_omega(2:end,:) - mesh_omega(1:end-1,:) ) ./ ...
 ( mesh_bigP(2:end,:) - mesh_bigP(1:end-1,:) );
contourf( bigXVals, cent(bigPVals), ddp_mesh_omega, 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "D/DP omega - rhoLin vs bigX, bigP" );

numFigs++; figure(numFigs);
contourf( bigXVals, bigPVals, sqrt(abs(mesh_omegaModel_rhoLin)), 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "omegaModel - rhoLin vs bigX, bigP" );

numFigs++; figure(numFigs);
ddx_mesh_omegaModel_rhoLin = ...
 ( mesh_omegaModel_rhoLin(:,2:end) - mesh_omegaModel_rhoLin(:,1:end-1) ) ./ ...
 ( mesh_bigX(:,2:end) - mesh_bigX(:,1:end-1) );
contourf( cent(bigXVals), bigPVals, ddx_mesh_omegaModel_rhoLin, 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "D/DX omegaModel - rhoLin vs bigX, bigP" );

numFigs++; figure(numFigs);
ddp_mesh_omegaModel_rhoLin = ...
 ( mesh_omegaModel_rhoLin(2:end,:) - mesh_omegaModel_rhoLin(1:end-1,:) ) ./ ...
 ( mesh_bigP(2:end,:) - mesh_bigP(1:end-1,:) );
contourf( bigXVals, cent(bigPVals), ddp_mesh_omegaModel_rhoLin, 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
plot( ...
  vecDeltaVals_rhoLin_levenberg(1,:)+bigX_guess, ...
  vecDeltaVals_rhoLin_levenberg(2,:)+bigP_guess, ...
  'wo-', 'markersize', 10 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "D/DP omegaModel - rhoLin vs bigX, bigP" );

numFigs++; figure(numFigs);
contourf( bigXVals, bigPVals, sqrt(abs(mesh_omegaModel_rhoSqQuad)), 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "omegaModel - rhoSqQuad vs bigX, bigP" );

numFigs++; figure(numFigs);
contourf( bigXVals, bigPVals, sqrt(abs(mesh_omegaModel_omegaQuad)), 21 );
colormap(mycmap(100));
hold on;
plot( bigX, bigP, 'w*', 'markersize', 25, 'linewidth', 3 );
plot( bigX_guess, bigP_guess, 'ws', 'markersize', 20, 'linewidth', 2 );
hold off;
grid on;
xlabel( "bigX" );
ylabel( "bigP" );
title( "omegaModel - omegaQuad vs bigX, bigP" );
