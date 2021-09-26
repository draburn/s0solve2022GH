function [ vecG, matH, retCode, datOut ] = findBestFit1D__calcLocalModel__tweak( dat_calcLocalModel, prm )
	%
	% Init
	commondefs;
	thisFile = "findBestFit1D__calcLocalModel__tweak";
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	omega = [];
	vecG = [];
	matH = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	%
	%
	% Using vecG via 2nd order finite differencing is usually good,
	% but not if we (somehow?) are near a local MAXimum.
	% In this case, we'll instead use one-sided finite differencing;
	% which side to use is somewhat arbitrary, but, we might as well
	% take the one which shows a greater reudction in omega.
	% Note that this necessarily corresponds to a large value of |dOmega/dZ|.
	omega = dat_calcLocalModel.omega;
	vecOmega_plus  = dat_calcLocalModel.vecOmega_plus;
	vecOmega_minus = dat_calcLocalModel.vecOmega_minus;
	vecEpsZ = dat_calcLocalModel.vecEpsZ;
	sizeZ = size(vecOmega_plus,1);
	gTweakCoeff = mygetfield( prm, "gTweakCoeff", 0.1 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega) );
		assert( 0.0 <= omega );
		assert( isrealarray(vecOmega_plus,[sizeZ,1]) );
		assert( isrealarray(vecOmega_minus,[sizeZ,1]) );
		assert( isrealarray(vecEpsZ,[sizeZ,1]) );
		assert( sum(vecEpsZ<=0.0)==0 );
		assert( isrealscalar(gTweakCoeff) );
		assert( 0.0 < gTweakCoeff );
		assert( gTweakCoeff <= 1.0 );
	end
	vecG_minus = (omega - vecOmega_minus)./(vecEpsZ);
	vecG_plus  = (vecOmega_plus  - omega)./(vecEpsZ);
	vecG_cent  = (vecOmega_plus  - vecOmega_minus)./(2.0*vecEpsZ);
	%
	vecG = zeros(sizeZ,1);
	for n=1:sizeZ
		%
		% DRaburn 2021.09.25:
		% I think this is right.
		% Hadn't properly tested before.
		% Hit for first time from test7_asym with setprngstates(6723760).
		%omega_plus  = vecOmega_plus(n);
		%omega_minus = vecOmega_minus(n);
		%
		plusIsGooderThan0  = false; % Unless...
		if ( ( abs(vecG_cent(n)) < gTweakCoeff*abs(vecG_plus(n)) ) ...
		  && ( omega_plus < omega ) )
			plusIsGooderThan0 = true;
		end
		%
		minusIsGooderThan0  = false; % Unless...
		if ( ( abs(vecG_cent(n)) < gTweakCoeff*abs(vecG_minus(n)) ) ...
		  && ( omega_minus < omega  ) )
			minusIsGooderThan0 = true;
		end
		%
		if ( plusIsGooderThan0 && minusIsGooderThan0 )
			if ( omega_plus <= omega_minus )
				whichToUse = +1;
			else
				whichToUse = -1;
			end
		elseif ( plusIsGooderThan0 )
			whichToUse = +1;
		elseif ( minusIsGooderThan0 )
			whichToUse = -1;
		else
			whichToUse = 0;
		end
		%
		if ( 0 == whichToUse )
			vecG(n) = vecG_cent(n);
		elseif ( -1 == whichToUse )
			% Apply gTweakCoeff* to prevent jump in value.
			vecG(n) = gTweakCoeff*vecG_minus(n);
		elseif ( +1 == whichToUse )
			% Apply gTweakCoeff* to prevent jump in value.
			vecG(n) = gTweakCoeff*vecG_plus(n);
		else
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "Invalid value of whichToUse (%g).", whichToUse ) );
			retCode = RETCODE__INTERNAL_INCONSISTENCY;
		end
	end
	%echo__vecG = vecG
	%
	%
	%
	% Find appropriate Hessian.
	% For starters, calc "H1" according to vecG from above...
	vecRho = dat_calcLocalModel.vecRho;
	matRho_plus  = dat_calcLocalModel.matRho_plus;
	matRho_minus = dat_calcLocalModel.matRho_minus;
	matOmegaDRho2 = dat_calcLocalModel.matOmegaDRho2;
	sizeRho = size(vecRho,1);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecRho,[sizeRho,1]) );
		assert( isrealarray(matRho_plus, [sizeZ,sizeRho]) );
		assert( isrealarray(matRho_minus,[sizeZ,sizeRho]) );
		assert( isrealarray(matOmegaDRho2,[sizeRho,sizeRho]) );
		assert( sum(sum((matOmegaDRho2'-matOmegaDRho2).^2)) <= eps150*sum(sum(matOmegaDRho2.^2)) );
	end
	%
	matRhoDZ = zeros(sizeZ,sizeRho);
	for n=1:sizeZ
		plusIsGooderThan0  = false; % Unless...
		if ( ( abs(vecG_cent(n)) < gTweakCoeff*abs(vecG_plus(n)) ) ...
		  && ( omega_plus < omega ) )
			plusIsGooderThan0 = true;
		end
		%
		minusIsGooderThan0  = false; % Unless...
		if ( ( abs(vecG_cent(n)) < gTweakCoeff*abs(vecG_minus(n)) ) ...
		  && ( omega_minus < omega  ) )
			minusIsGooderThan0 = true;
		end
		%
		if ( plusIsGooderThan0 && minusIsGooderThan0 )
			if ( omega_plus <= omega_minus )
				whichToUse = +1;
			else
				whichToUse = -1;
			end
		elseif ( plusIsGooderThan0 )
			whichToUse = +1;
		elseif ( minusIsGooderThan0 )
			whichToUse = -1;
		else
			whichToUse = 0;
		end
		%
		if ( 0 == whichToUse )
			matRhoDZ(n,:) = ( matRho_plus(n,:) - matRho_minus(n,:) ) / ( 2.0*vecEpsZ(n) );
		elseif ( -1 == whichToUse )
			matRhoDZ(n,:) = ( vecRho' - matRho_minus(n,:) ) / vecEpsZ(n);
		elseif ( +1 == whichToUse )
			matRhoDZ(n,:) = ( matRho_plus(n,:) - vecRho' ) / vecEpsZ(n);
		else
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "Invalid value of whichToUse (%g).", whichToUse ) );
			retCode = RETCODE__INTERNAL_INCONSISTENCY;
		end
	end
	%
	matH = matRhoDZ * matOmegaDRho2 * (matRhoDZ');
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(matH,[sizeZ,sizeZ]) );
		assert( sum(sum((matH'-matH).^2)) <= eps150*sum(sum(matH.^2)) );
	end
	%
	% For regularization, we could consider both "H1" and "H2".
	% Ultimately, however, we'd need to use some decisive regularization anyway,
	% so, we might as well just do a decisive algorithm.
	% First, force all elements on main diagonal to be clearly positive.
	hScale = sqrt(sum( sum(matH.^2)/sizeZ ));
	if ( 0.0 == hScale )
		% Well,... damn.
		% Take matH to be the (pos-def) diag matrix s.t. (model) omega min is zero...
		if ( valLev >= VALLEV__MEDIUM )
		if ( 0.0 == omega )
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		end
		mu = sqrt(sum(vecG.^2))/(2.0*omega);
		matH = mu*eye(sizeZ,sizeZ);
		%
		retCode = RETCODE__SUCCESS;
		return;
	end
	for n=1:sizeZ
	if ( matH(n,n) <= eps075*hScale )
		matH(n,n) = eps075*hScale;
	end
	end
	%
	%
	%
	% Now, make sure matH is invertible / has a valid Cholesky factorization.
	mu0 = mygetfield( prm, "mu0", eps050 );
	mu1 = mygetfield( prm, "mu1", 1.0./eps050 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(mu0) );
		assert( isrealscalar(mu1) );
		assert( isrealscalar(muStep) );
		assert( 0.0 < mu0 );
		assert( mu0 < mu1 );
		assert( 1.0 < muStep );
	end
	%
	iterCount = 0;
	errFlag = true;
	while (errFlag)
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
			  "Failed to regularize Hessian (%g).", mu_trial ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		%
		matH_trial = matH + mu_trial*hScale*eye(sizeZ,sizeZ);
		[ foo, errFlag ] = chol(matH_trial);
	end
	matH = matH_trial;
	%echo__matH = matH
	%
	if (false)
	% Make sure that the minimum omega according to the model is >= zero?
	hScale = sqrt(sum( sum(matH.^2)/sizeZ )); % Updated!
	iterCount = 0;
	errFlag = true;
	while (errFlag)
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
			  "Failed to make omega model well-behaved (%g).", mu_trial ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		%
		matH_trial = matH + mu_trial*hScale*eye(sizeZ,sizeZ);
		vecDelta_trial = -matH_trial \ vecG;
		omega_trial = omega + vecDelta_trial'*vecG + 0.5*vecDelta_trial'*matH_trial*vecDelta_trial;
		if ( omega_trial < -eps075*abs(omega) )
			errFlag = true;
		else
			errFlag = false;
		end
	end
	matH = matH_trial;
	end
	%
	retCode = RETCODE__SUCCESS;
	return;
end
