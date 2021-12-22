function [ vecDelta, retCode, datOut ] = findBestFit1D__findStep( funchRho, rhoArgs, vecZ, prm )
	%
	% Init
	commondefs;
	thisFile = "findBestFit1D__findStep";
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	vecDelta = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	%
	%
	sizeZ = max(size(vecZ));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecZ,[sizeZ,1]) );
	end
	%
	[ omega, vecG, matH, retCode, datOut_calcLocalModel ] = findBestFit1D__calcLocalModel( ...
	  funchRho, rhoArgs, vecZ, prm );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcLocalModel() returned %s.", retcode2str(retCode) ) );
		% Leave retCode as is.
		return;
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega) );
		assert( 0.0 <= omega );
		assert( isrealarray(vecG,[sizeZ,1]) );
		assert( isrealarray(matH,[sizeZ,sizeZ]) );
	end
	datOut.calcLocalModel = datOut_calcLocalModel;
	vecRho = datOut_calcLocalModel.vecRho;
	sizeRho = size(vecRho,1);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealvector(vecRho,[sizeRho,1]) );
	end
	%
	% This should be passed in...
	if ( mygetfield(prm,"useCustomOmega",false) )
		funchOmega      = mygetfield(prm,"funchOmega");
		funchOmegaDRho  = mygetfield(prm,"funchOmegaDRho");
		funchOmegaDRho2 = mygetfield(prm,"funchOmegaDRho2");
	else
		funchOmega      = @(rho)( 0.5 * sum(rho.^2) );
		funchOmegaDRho  = @(rho)( rho );
		funchOmegaDRho2 = @(rho)( eye(sizeRho,sizeRho) );
	end
	%
	%
	%
	% Backtrack until we find an acceptable delta or give up.
	mu0 = mygetfield( prm, "mu0", eps050 );
	mu1 = mygetfield( prm, "mu1", 1.0./eps050 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	requireOmegaModelToBePositive = mygetfield( prm, "requireOmegaModelToBePositive", false );
	sufficientDecreaseCoeff = mygetfield( prm, "sufficientDecreaseCoeff", 0.01 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
		assert( issize(requireOmegaModelToBePositive,[1,1]) );
		assert( isrealscalar(sufficientDecreaseCoeff) );
		assert( 0.0 < sufficientDecreaseCoeff );
		assert( sufficientDecreaseCoeff < 1.0 );
	end
	%
	iterCount = 0;
	while (true)
		retCode = RETCODE__NOT_SET;
		iterCount++;
		switch (iterCount)
		case 1
			mu_trial = 0.0;
		case 2
			mu_trial = mu0;
		otherwise
			mu_trial *= muStep;
		end
		if ( mu_trial > mu1 )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to find a good step (%g).", mu_trial ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		[ vecDelta_trial, retCode, datOut_calcDelta_trial ] = findBestFit1D__findStep__calcDelta( ...
		  vecZ, omega, vecG, matH, mu_trial, prm );
		if ( RETCODE__SUCCESS ~= retCode )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "__calcDelta() returned %s.", retcode2str(retCode) ) );
			continue;
		end
		retCode = RETCODE__NOT_SET;
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealarray(vecDelta_trial,[sizeZ,1]) );
		end
		%
		omegaModel_trial = omega + vecDelta_trial' * vecG + 0.5 * vecDelta_trial' * matH * vecDelta_trial;
		if ( omegaModel_trial > omega )
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "Generated delta increases omega model (%g).", mu_trial ) );
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "Omegas: %g - %g = %g.", omegaModel_trial, omega, omegaModel_trial - omega ) );
			retCode = RETCODE__INTERNAL_INCONSISTENCY;
			return;
		end
		if ( requireOmegaModelToBePositive )
		if ( omegaModel_trial < 0.0 )
			continue;
		end
		end
		%
		%
		%
		[ errFlag, vecRho_trial ] = funchRho( rhoArgs, vecZ + vecDelta_trial );
		if (errFlag)
			continue;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealarray(vecRho_trial,[sizeRho,1]) );
		end
		%
		omega_trial = funchOmega( vecRho_trial );
		if (~isrealscalar(omega_trial))
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchOmega() failed for trial (not a real scalar) (%g).", mu_trial ) );
			continue;
		end
		if ( omega_trial < 0.0 )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchOmega() failed for trial (negative) (%g).", mu_trial ) );
			continue;
		end
		%
		if ( omega_trial > omega )
			continue;
		end
		if ( abs(omega_trial - omega) < sufficientDecreaseCoeff*abs(omegaModel_trial - omega) );
			continue;
		end
		%
		%
		% Success!
		vecDelta = vecDelta_trial;
		retCode = RETCODE__SUCCESS;
		return;
	end
end
