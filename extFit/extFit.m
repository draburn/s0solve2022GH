function [ datOut, retCode ] = extFit( bigX0, bigP0, rvecX, rvecF, rvecW=[], prm=[] )
	%
	commondefs; thisFile = "extFit";
	msg( thisFile, __LINE__, "Note that I may be better than diag(diag(H))." );
	msg( thisFile, __LINE__, "TODO: Include bounds in X!" );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	datOut = [];
	retCode = RETCODE__NOT_SET;
	%
	if (isempty(rvecW))
		rvecW = ones(size(rvecX));
	end
	numPts = size(rvecX,2);
	assert( isrealscalar(bigX0) );
	assert( isrealscalar(bigP0) );
	assert( isrealarray(rvecX,[1,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(rvecW,[1,numPts]) );
	%
	datOut.bigX0 = bigX0;
	datOut.bigP0 = bigP0;
	datOut.rvecX = rvecX;
	datOut.rvecF = rvecF;
	datOut.rvecW = rvecW;
	datOut.prm = prm;
	%
	%dat0 = extFit_calcOmega( bigX0, bigP0, rvecX, rvecF, rvecW );
	%omega0 = dat0.omega0;
	%
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	assert( isrealscalar(iterLimit) );
	assert( iterLimit >= 1 );
	%
	omegaTol = mygetfield( prm, "omegaTol", 0.5*eps*sum(rvecF.^2) );
	assert( isrealscalar(omegaTol) );
	assert( omegaTol > 0.0 );
	%
	deltaXTol = sqrt(eps)*(max(rvecX)-min(rvecX));
	deltaXTol = mygetfield( prm, "deltaXTol", deltaXTol );
	assert( isrealscalar(deltaXTol) );
	assert( deltaXTol > 0.0 );
	deltaPTol = mygetfield( prm, "deltaPTol", sqrt(sqrt(eps)) );
	assert( isrealscalar(deltaPTol) );
	assert( deltaPTol > 0.0 );
	%
	bigX = bigX0;
	bigP = bigP0;
	iterCount = 0;
	%
	prm_calcGradHess = mygetfield( prm, "prm_calcGradHess", [] );
	dat_calcGradHess = extFit_calcGradHess( bigX, bigP, rvecX, rvecF, rvecW, prm_calcGradHess );
	omega = dat_calcGradHess.omega0;
	if ( omega <= omegaTol )
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Initial guess is converged (%g <= %g).", omega, omegaTol ) );
		retCode = RETCODE__SUCCESS;
		datOut.bigX = bigX;
		datOut.bigP = bigP;
		return;
	end
	vecG = dat_calcGradHess.vecG;
	matH = dat_calcGradHess.matH1;
	if ( sum(diag(matH)<=0.0) != 0 )
		msg_warn( verbLev, thisFile, __LINE__, "Hessian diagonal has a non-positive element." );
		retCode = RETCODE__ALGORITHM_BREAKDOWN;
		datOut.bigX = bigX;
		datOut.bigP = bigP;
		return;
	end
	matD = diag(abs(diag(matH)));
	%matD = eye(2,2);
	mu_trial = 0.0;
	while (1)
		iterCount++;
		trialIsValid = true;
		deltaVec_trial = - ( matH + mu_trial*matD )	\ vecG;
		bigX_trial = bigX + deltaVec_trial(1);
		bigP_trial = bigP + deltaVec_trial(2);
		if ( bigP_trial <= 0.0 )
			trialIsValid = false;
		end
		if ( trialIsValid )
			omega_trial = extFit_calcOmega( rvecX, rvecF, bigX_trial, bigP_trial, rvecW );
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
				  " %3d;  %11.3e, %11.3e, %11.3e;  %11.3e, %11.3e, %11.3e; %2d, %11.3e, %11.3e, %11.3e; %11.3e", ...
			  iterCount, ...
			  bigX, ...
			  bigP, ...
			  omega, ...
			  mu_trial, ...
			  deltaVec_trial(1), ...
			  deltaVec_trial(2), ...
			  trialIsValid, ...
			  bigX_trial, ...
			  bigP_trial, ...
			  omega_trial, ...
			  omega_trial - omega ) );
		else
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
				  "%3d;  %11.3e, %11.3e, %11.3e;  %11.3e, %11.3e, %11.3e; %2d, %11.3e, %11.3e, %11.3e; %11.3e", ...
			  iterCount, ...
			  bigX, ...
			  bigP, ...
			  omega, ...
			  mu_trial, ...
			  deltaVec_trial(1), ...
			  deltaVec_trial(2), ...
			  trialIsValid, ...
			  bigX_trial, ...
			  bigP_trial, ...
			  -1.0, ...
			  -1.0 ) );
		end
		%
		if ( trialIsValid && omega_trial <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged (%g <= %g).", omega_trial, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			datOut.bigX = bigX;
			datOut.bigP = bigP;
			return;
		elseif ( abs(deltaVec_trial(1)) <= deltaXTol && abs(deltaVec_trial(2)) <= deltaPTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Found local extermum (|%g| <= %g, |%g| <= %g).", ...
			  deltaVec_trial(1), deltaXTol, deltaVec_trial(2), deltaPTol ) );
			retCode = RETCODE__SUCCESS;
			datOut.bigX = bigX;
			datOut.bigP = bigP;
			return;
		elseif ( trialIsValid && omega_trial < omega )
			% Move, update grad & hess, and continue.
			bigX = bigX_trial;
			bigP = bigP_trial;
			datOut.bigX = bigX;
			datOut.bigP = bigP;
			%
			dat_calcGradHess = extFit_calcGradHess( bigX, bigP, rvecX, rvecF, rvecW, prm_calcGradHess );
			omega = dat_calcGradHess.omega0;
			vecG = dat_calcGradHess.vecG;
			matH = dat_calcGradHess.matH1;
			if ( sum(diag(matH)<=0.0) != 0 )
				msg_warn( verbLev, thisFile, __LINE__, "Hessian diagonal has a non-positive element." );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				datOut.bigX = bigX;
				datOut.bigP = bigP;
				return;
			end
			matD = diag(abs(diag(matH)));
			%matD = eye(2,2);
			mu_trial = 0.0;
			continue;
		elseif ( iterCount >= iterLimit )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached iterLimit (%d).", iterLimit ) );
			retCode = RETCODE__IMPOSED_STOP;
			datOut.bigX = bigX;
			datOut.bigP = bigP;
			return;
		else
			% Increase mu and try again.
			if ( 0.0 == mu_trial )
				mu_trial = 0.1;
			else
				mu_trial *= 10.0;
			end
			continue;
		end
	end
	%
end
