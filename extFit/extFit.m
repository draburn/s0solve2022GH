function [ datOut, retCode ] = extFit( bigX0, bigP0, rvecX, rvecF, rvecW=[], prm=[] )
	clear;
	%setprngstates();
	%setprngstates(19719664); % eye(2,2) may be better here.
	%setprngstates(77173824); % Better initial guess would help; goes negative in P.
	setprngstates(80489888);
	numPts = 5 + round(abs(randn()*exp(abs(randn()))))
	bigX_secret = randn()*exp(abs(3.0*randn()))
	bigP_secret = 1.0 + 3.0*abs(randn())
	bigA_secret = randn()*exp(abs(3.0*randn()));
	bigB_secret = randn()*exp(abs(3.0*randn()));
	rvecX = sort([ ...
	  bigX_secret - abs(randn(1,2)), ...
	  bigX_secret + abs(randn(1,2)), ...
	  bigX_secret + randn(1,numPts-4) ]);
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ).^bigP_secret );
	rvecF = funchF(rvecX);
	rvecW = [];
	prm = [];
	index0 = 1;
	if ( bigB_secret > 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) > rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	elseif ( bigB_secret < 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) < rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	end
	bigX0 = (rvecX(index0+1)+rvecX(index0-1))/2.0
	bigP0 = 2.0
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	commondefs; thisFile = "extFit";
	msg( thisFile, __LINE__, "Note that I may be better than diag(diag(H))." );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
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
	omegaTol = mygetfield( prm, "omegaTol", 0.5*eps*sum(rvecF.^2) )
	assert( isrealscalar(omegaTol) );
	assert( omegaTol > 0.0 );
	%
	deltaXTol = sqrt(eps)*(max(rvecX)-min(rvecX));
	deltaXTol = mygetfield( prm, "deltaXTol", deltaXTol )
	assert( isrealscalar(deltaXTol) );
	assert( deltaXTol > 0.0 );
	deltaPTol = mygetfield( prm, "deltaPTol", sqrt(sqrt(eps)) )
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
			dat_trial = extFit_calcOmega( bigX_trial, bigP_trial, rvecX, rvecF, rvecW );
			omega_trial = dat_trial.omega;
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
