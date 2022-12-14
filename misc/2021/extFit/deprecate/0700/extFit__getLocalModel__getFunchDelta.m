function datOut = extFit__getLocalModel__getFunchDelta( vecG, matH, prm = [] )
	thisFile = "extFit__getLocalModel__getFunchDelta";
	msg( thisFile, __LINE__, "DEPRECATED." );
	% We should maybe do something if matH is not pos-def?
	%
	matI = eye(2,2);
	vecDeltaNewt = -matH \ vecG;
	matS = diag(abs(diag(matH)));
	matHMS = matH - matS;
	matHMI = matH - matI;
	vecGS = matS \ vecG;
	vecGHat = -vecG * norm(vecDeltaNewt)/norm(vecG);
	vecGSHat = -vecGS * norm(vecDeltaNewt)/norm(vecGS);
	%
	datOut.newton = @(lambda)( lambda*vecDeltaNewt );
	datOut.gradDir = @(lambda)( lambda*vecGHat );
	datOut.gradSclDir = @(lambda)( lambda*vecGSHat );
	datOut.levenberg = @(lambda)( -lambda*( (lambda*matHMI+matI)\vecG ) );
	datOut.levMarq = @(lambda)( -lambda*( (lambda*matHMS+matS)\vecG ) );
	% gradCurve?
	% gradScaledCurve?
return;
end

%!test
%!	thisFile = "test extFit__getLocalModel__getFunchDelta";
%!	bigA = 1.0;
%!	bigB = 1.0;
%!	bigX = 0.5;
%!	bigP = 2.0;
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = linspace( -2.0, 3.0, 6 );
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX + 0.1;
%!	bigP_guess = bigP + 0.1;
%!	
%!	dat_localModel_rhoLin    = extFit__getLocalModel_rhoLin(    bigX_guess, bigP_guess, xVals, fVals );
%!	dat_localModel_rhoSqQuad = extFit__getLocalModel_rhoSqQuad( bigX_guess, bigP_guess, xVals, fVals );
%!	dat_localModel_omegaQuad = extFit__getLocalModel_omegaQuad( bigX_guess, bigP_guess, xVals, fVals );
%!	lambdaVals = linspace(0.0,1.0,101);
%!	n = 0;
%!	for lambda=lambdaVals
%!		n++;
%!		%
%!		vecDeltaVals_rhoLin_newton(:,n)     = dat_localModel_rhoLin.dat_funchDelta.newton(     lambda );
%!		vecDeltaVals_rhoLin_gradDir(:,n)    = dat_localModel_rhoLin.dat_funchDelta.gradDir(    lambda );
%!		vecDeltaVals_rhoLin_gradSclDir(:,n) = dat_localModel_rhoLin.dat_funchDelta.gradSclDir( lambda );
%!		vecDeltaVals_rhoLin_levenberg(:,n)  = dat_localModel_rhoLin.dat_funchDelta.levenberg(  lambda );
%!		vecDeltaVals_rhoLin_levMarq(:,n)    = dat_localModel_rhoLin.dat_funchDelta.levMarq(    lambda );
%!		%
%!		vecDeltaVals_rhoSqQuad_newton(:,n)     = dat_localModel_rhoSqQuad.dat_funchDelta.newton(     lambda );
%!		vecDeltaVals_rhoSqQuad_gradDir(:,n)    = dat_localModel_rhoSqQuad.dat_funchDelta.gradDir(    lambda );
%!		vecDeltaVals_rhoSqQuad_gradSclDir(:,n) = dat_localModel_rhoSqQuad.dat_funchDelta.gradSclDir( lambda );
%!		vecDeltaVals_rhoSqQuad_levenberg(:,n)  = dat_localModel_rhoSqQuad.dat_funchDelta.levenberg(  lambda );
%!		vecDeltaVals_rhoSqQuad_levMarq(:,n)    = dat_localModel_rhoSqQuad.dat_funchDelta.levMarq(    lambda );
%!		%
%!		vecDeltaVals_omegaQuad_newton(:,n)     = dat_localModel_omegaQuad.dat_funchDelta.newton(     lambda );
%!		vecDeltaVals_omegaQuad_gradDir(:,n)    = dat_localModel_omegaQuad.dat_funchDelta.gradDir(    lambda );
%!		vecDeltaVals_omegaQuad_gradSclDir(:,n) = dat_localModel_omegaQuad.dat_funchDelta.gradSclDir( lambda );
%!		vecDeltaVals_omegaQuad_levenberg(:,n)  = dat_localModel_omegaQuad.dat_funchDelta.levenberg(  lambda );
%!		vecDeltaVals_omegaQuad_levMarq(:,n)    = dat_localModel_omegaQuad.dat_funchDelta.levMarq(    lambda );
%!	end
%!	
%!	numFigs = 0;
%!	
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  lambdaVals, vecDeltaVals_rhoLin_newton(1,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_rhoLin_levenberg(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_rhoLin_levMarq(1,:), '^-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_newton(1,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(1,:), '^-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_newton(1,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_levenberg(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_levMarq(1,:), '^-' );
%!	grid on;
%!	xlabel( "lambda" );
%!	ylabel( "deltaX" );
%!	legend( ...
%!	  "rhoLin - newt", ...
%!	  "rhoLin - lev", ...
%!	  "rhoLin - levMarq", ...
%!	  "rhoSqQuad - newt", ...
%!	  "rhoSqQuad - lev", ...
%!	  "rhoSqQuad - levMarq", ...
%!	  "omegaQuad - newt", ...
%!	  "omegaQuad - lev", ...
%!	  "omegaQuad - levMarq", ...
%!	  "location", "southwestoutside" );
%!	
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  lambdaVals, vecDeltaVals_rhoLin_newton(2,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_rhoLin_levenberg(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_rhoLin_levMarq(2,:), '^-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_newton(2,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_levenberg(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_rhoSqQuad_levMarq(2,:), '^-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_newton(2,:), '*-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_levenberg(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_omegaQuad_levMarq(2,:), '^-' );
%!	grid on;
%!	xlabel( "lambda" );
%!	ylabel( "deltaP" );
%!	legend( ...
%!	  "rhoLin - newt", ...
%!	  "rhoLin - lev", ...
%!	  "rhoLin - levMarq", ...
%!	  "rhoSqQuad - newt", ...
%!	  "rhoSqQuad - lev", ...
%!	  "rhoSqQuad - levMarq", ...
%!	  "omegaQuad - newt", ...
%!	  "omegaQuad - lev", ...
%!	  "omegaQuad - levMarq", ...
%!	  "location", "southwestoutside" );
