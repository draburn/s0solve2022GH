function [ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm=[] );
	%
	% Do standard stuff.
	commondefs;
	thisFile = "findBestFit1D";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	vecZ = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	%
	% Unpack & validate primary data.
	sizeZ = max(size(vecZ0));
	if ( valLev >= VALLEV__LOW )
		assert( isrealvector(vecZ0,sizeZ) );
	end
	%
	[ errFlag, vecRho0 ] = funchRho( rhoArgs, vecZ0 );
	if ( errFlag )
		msg_error( verbLev, thisFile, __LINE__, "funchRho() failed for initial guess." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	sizeRho = max(size(vecRho0));
	if ( valLev >= VALLEV__LOW )
		assert( isrealvector(vecRho0,sizeRho) );
	end
	%
	if ( mygetfield(prm,"useCustomOmega",false) )
		funchOmega = mygetfield(prm,"funchOmega");
	else
		funchOmega = @(rho)( 0.5 * sum(rho.^2) );
	end
	omega0 = funchOmega(vecRho0);
	if ( valLev >= VALLEV__LOW )
		assert( isrealscalar(omega0) );
		assert( omega0 >= 0.0 );
	end
	%
	%
	% Main loop.
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	omegaTol = mygetfield( prm, "omeagTol", eps050 );
	deltaOmegaRelTol = mygetfield( prm, "deltaOmegaRelTol", eps025 );
	vecDeltaTol = mygetfield( prm, "vecDeltaTol", eps025*ones(size(vecZ0)) );
	%
	vecZ = vecZ0;
	omega = omega0;
	iterCount = 0;
	msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
	  "   %3d;   %11s,  %11.3e;   %11s,  %11.3e.", iterCount, "-1", 0.0, "-1", omega ) );
	while (1)
		%
		%
		% Check pre-trial stopping criteria.
		if ( ~isempty(omegaTol) )
		if ( omegaTol >= 0.0 )
		if ( omega <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached omegaTol (%0.3e <= %0.3e).", omega, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			return;
		end
		end
		end
		%
		if ( ~isempty(iterLimit) )
		if ( iterLimit >= 0 )
		if ( iterCount >= iterLimit )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached iterLimit ( %d >= %d ).", iterCount, iterLimit ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		end
		end
		%
		%
		%
		% Find next trial.
		[ vecDelta, retCode, datOut_findStep ] = findBestFit1D__findStep( funchRho, rhoArgs, vecZ, prm );
		if ( RETCODE__SUCCESS != retCode )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "__findStep() returned %s.", retcode2str(retCode) ) );
			% Leave retCode unmodified (for now).
			return;
		end
		retCode = RETCODE__NOT_SET;
		if ( valLev >= VALLEV__MEDIUM )
			%assert( isrealvector(vecDelta,sizeZ) ); % Not good enough.
			assert( isrealarray(vecDelta,size(vecZ)) );
		end
		%
		vecZ_trial = vecZ + vecDelta;
		[ errFlag_trial, vecRho_trial ] = funchRho( rhoArgs, vecZ_trial );
		if ( errFlag_trial )
			msg_error( verbLev, thisFile, __LINE__, "funchRho() failed for step generated by __findStep." );
			retCode = RETCODE__UNSPECIFIC_ERROR;
			return;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealvector(vecRho_trial,sizeRho) );
		end
		%
		omega_trial = funchOmega( vecRho_trial );
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealscalar(omega_trial) );
			assert( 0.0 <= omega_trial );
		end
		%
		if ( omega_trial >= omega )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega ( %0.3e - %0.3e = %0.3e ).", omega_trial, omega, omega_trial-omega ) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return;
		end
		%
		%
		%
		% Check post-trial stopping criteria.
		if ( ~isempty(deltaOmegaRelTol) )
		if ( deltaOmegaRelTol >= 0.0 )
		if ( abs(omega-omega_trial) < deltaOmegaRelTol*abs(omega) )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to decrease omega sufficiently ( %0.3e - %0.3e = %0.3e ).", omega_trial, omega, omega_trial-omega ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		end
		end
		%
		if ( ~isempty(vecDeltaTol) )
			amOverDeltaTol = false;
			for n=1:sizeZ
			if ( vecDeltaTol(n) >= 0.0 )
			if ( abs(vecZ_trial(n)-vecZ(n)) > vecDeltaTol(n) )
				amOverDeltaTol = true;
				break;
			end
			end
			end
			if (~amOverDeltaTol)
				msg_main( verbLev, thisFile, __LINE__, sprintf( ...
				  "Reached deltaTol in every component ( %0.3e ).", norm(vecDelta) ) );
				retCode = RETCODE__SUCCESS;
				return;
			end
		end
		%
		%
		%
		% Iterate.
		iterCount++;
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
		  "  %3d;   %11.3e,  %11.3e;   %11.3e,  %11.3e.", ...
		  iterCount, norm(vecDelta), norm(vecZ_trial-vecZ0), omega-omega_trial, omega_trial ) );
		vecZ = vecZ_trial;
		vecRho = vecRho_trial;
		omega = omega_trial;
	end
end


%!test
%!	commondefs;
%!	setprngstates(0);
%!	f0 = 0.2;
%!	f1 = 0.4;
%!	x0 = 0.6;
%!	p = 3.5;
%!	rhoArgs.xVals = linspace( -1.0, 2.0, 5 );
%!	rhoArgs.fVals = f0+f1*abs(rhoArgs.xVals-x0).^p;
%!	rhoArgs.dVals = ones(size(rhoArgs.xVals));
%!	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
%!	%
%!	vecZ = [ 0.0, 2.0 ];
%!	[ errFlag, vecRho ] = funchRho( rhoArgs, vecZ );
%!	echo__errFlag = errFlag;
%!	echo__vecRho = vecRho;
%!	%
%!	vecZ = [ x0, p ];
%!	[ errFlag, vecRho ] = funchRho( rhoArgs, vecZ );
%!	echo__errFlag = errFlag;
%!	echo__vecRho = vecRho;
%!	%
%!	prm = [];
%!	echo__prm = prm
%!	vecZ0 = [ 0.0, 2.0 ]
%!	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm=[] );
%!	echo__vecZ = vecZ
%!	echo__retCode = retCode;
%!	echo__datOut = datOut;
%!	assert( retCode == RETCODE__SUCCESS );

%!test
%!	commondefs;
%!	setprngstates(0);
%!	%
%!	rhoArgs.xVals = linspace( 0.0, 2.0, 3 );
%!	rhoArgs.fVals = abs(rhoArgs.xVals).^3;
%!	rhoArgs.dVals = ones(size(rhoArgs.xVals));
%!	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
%!	%
%!	prm = [];
%!	echo__prm = prm
%!	vecZ0 = [ 0.0, 2.0 ]
%!	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm=[] );
%!	echo__vecZ = vecZ
%!	echo__retCode = retCode;
%!	echo__datOut = datOut;
%!	assert( retCode == RETCODE__SUCCESS );
