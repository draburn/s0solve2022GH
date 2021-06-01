thisFile = "test_extFit__getLocalModel";
bigA = 1.0;
bigB = 1.0;
bigX = 0.5;
bigP = 2.0;
funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
xVals = linspace( -2.0, 3.0, 6 );
fVals = funchF(xVals);
bigX_guess = bigX + 0.0;
bigP_guess = bigP + 0.01;
%
dat_localModel_rhoLin    = extFit__getLocalModel_rhoLin(    bigX_guess, bigP_guess, xVals, fVals );
dat_localModel_rhoSqQuad = extFit__getLocalModel_rhoSqQuad( bigX_guess, bigP_guess, xVals, fVals );
dat_localModel_omegaQuad = extFit__getLocalModel_omegaQuad( bigX_guess, bigP_guess, xVals, fVals );
lambdaVals = linspace(0.0,1.0,101);
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
end

numFigs = 0;

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, vecDeltaVals_rhoLin_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(1,:), '^-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(1,:), '^-', ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(1,:), '^-' );
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
  lambdaVals, vecDeltaVals_rhoLin_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(1,:), '^-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(1,:), '^-', ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(1,:), '*-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(1,:), 'o-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(1,:), '^-' );
grid on;
xlabel( "lambda" );
ylabel( "deltaX" );

numFigs++; figure(numFigs);
plot( ...
  lambdaVals, vecDeltaVals_rhoLin_newton(2,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoLin_levenberg(2,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoLin_levMarq(2,:), '^-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_newton(2,:), '*-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(2,:), 'o-', ...
  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(2,:), '^-', ...
  lambdaVals, vecDeltaVals_omegaQuad_newton(2,:), '*-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levenberg(2,:), 'o-', ...
  lambdaVals, vecDeltaVals_omegaQuad_levMarq(2,:), '^-' );
grid on;
xlabel( "lambda" );
ylabel( "deltaP" );
