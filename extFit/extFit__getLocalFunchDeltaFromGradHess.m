function datOut = extFit__getLocalFunchDeltaFromGradHess( vecG, matH, prm = [] )
	thisFile = "extFit__getLocalFunchDeltaFromGradHess";
	% We should maybe do something if matH is not pos-def?
	%
	matI = eye(2,2);
	vecDeltaNewt = matH \ vecG;
	matS = diag(abs(diag(matH)));
	matHMS = matH - matS;
	matHMI = matH - matI;
	vecGS = matS \ vecG;
	vecGHat = vecG * norm(vecDeltaNewt)/norm(vecG);
	vecGSHat = vecGS * norm(vecDeltaNewt)/norm(vecGS);
	%
	datOut.funch_delta_newt = @(lambda)( lambda*vecDeltaNewt );
	datOut.funch_delta_gradDir = @(lambda)( lambda*vecGHat );
	datOut.funch_delta_gradSclDir = @(lambda)( lambda*vecGSHat );
	datOut.funch_delta_lev = @(lambda)( lambda*( (lambda*matHMI+matI)\vecG ) );
	datOut.funch_delta_levMarq = @(lambda)( lambda*( (lambda*matHMS+matS)\vecG ) );
	% gradCurve?
	% gradScaledCurve?
return;
end

%!test
%!	thisFile = "test extFit__getFunchDelta"
%!	bigA = 1.0;
%!	bigB = 1.0;
%!	bigX = 0.5;
%!	bigP = 2.0;
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = linspace( -2.0, 3.0, 6 );
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX + 0.1;
%!	bigP_guess = bigP + 0.1;
%!	dat_locGH_rhoLin = extFit__getLocalGradHess_rhoLin( bigX_guess, bigP_guess, xVals, fVals );
%!	dat_locFunchDelta_rhoLin = extFit__getLocalFunchDeltaFromGradHess( dat_locGH_rhoLin.vecG, dat_locGH_rhoLin.matH );
%!	lambdaVals = linspace(0.0,1.0,101);
%!	n = 0;
%!	for lambda=lambdaVals
%!		n++;
%!		vecDeltaVals_newt_rhoLin(:,n) = dat_locFunchDelta_rhoLin.funch_delta_newt( lambda );
%!		vecDeltaVals_gradDir_rhoLin(:,n) = dat_locFunchDelta_rhoLin.funch_delta_gradDir( lambda );
%!		vecDeltaVals_gradSclDir_rhoLin(:,n) = dat_locFunchDelta_rhoLin.funch_delta_gradSclDir( lambda );
%!		vecDeltaVals_lev_rhoLin(:,n) = dat_locFunchDelta_rhoLin.funch_delta_lev( lambda );
%!		vecDeltaVals_levMarq_rhoLin(:,n) = dat_locFunchDelta_rhoLin.funch_delta_levMarq( lambda );
%!	end
%!
%!	numFigs = 0;
%!
%!	numFigs++; figure(numFigs);	
%!	plot( ...
%!	  lambdaVals, vecDeltaVals_newt_rhoLin(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_gradDir_rhoLin(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_gradSclDir_rhoLin(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_lev_rhoLin(1,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_levMarq_rhoLin(1,:), 'o-' );
%!	grid on;
%!	xlabel( "lambda" );
%!	ylabel( "deltaX" );
%!	legend( ...
%!	  "newt - rhoLin", ...
%!	  "gradDir - rhoLin", ...
%!	  "gradSclDir - rhoLin", ...
%!	  "lev - rhoLin", ...
%!	  "levMarq - rhoLin", ...
%!	  "location", "southwest" );
%!
%!	numFigs++; figure(numFigs);	
%!	plot( ...
%!	  lambdaVals, vecDeltaVals_newt_rhoLin(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_gradDir_rhoLin(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_gradSclDir_rhoLin(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_lev_rhoLin(2,:), 'o-', ...
%!	  lambdaVals, vecDeltaVals_levMarq_rhoLin(2,:), 'o-' );
%!	grid on;
%!	xlabel( "lambda" );
%!	ylabel( "deltaP" );
%!	legend( ...
%!	  "newt - rhoLin", ...
%!	  "gradDir - rhoLin", ...
%!	  "gradSclDir - rhoLin", ...
%!	  "lev - rhoLin", ...
%!	  "levMarq - rhoLin", ...
%!	  "location", "southwest" );
%!
