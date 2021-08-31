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
	
	% DRaburn 2021.08.29.
	% Dev hack!
	[ foo1, foo2, foo3, retCode, datOut_calcLocalModel ] = findBestFit1D__calcLocalModel( ...
	  funchRho, rhoArgs, vecZ, prm );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcLocalModel() returned %s.", retcode2str(retCode) ) );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	clear foo1;
	clear foo2;
	clear foo3;
	[ omega, vecG, matH, retCode, datOut_tweakLocalModel ] = findBestFit1D__tweakLocalModel( ...
	  datOut_calcLocalModel, prm );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__tweakLocalModel() returned %s.", retcode2str(retCode) ) );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	
	%
	%
	sizeZ = max(size(vecZ));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealvector(vecZ,sizeZ) );
	end
	%
	%
	%
if (0)
	% Calculate Hessian.
	[ errFlag, vecRho0 ] = funchRho( rhoArgs, vecZ );
	if (errFlag)
		msg_error( verbLev, thisFile, __LINE__, "funchRho() failed for current guess." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	sizeRho = max(size(vecRho0));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealvector(vecRho0,sizeRho) );
	end
	%
	vecEpsZ = mygetfield( prm, "vecEpsZ", eps050*ones(size(vecZ)) );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealvector(vecEpsZ,sizeZ) );
		for n=1:sizeZ
			assert( vecEpsZ(n) > 0.0 );
		end
	end
	%
	matRhoP = zeros(sizeRho,sizeZ);
	for n=1:sizeZ
		epsZ_temp = vecEpsZ(n);
		%
		vecZ_temp = vecZ;
		vecZ_temp(n) = vecZ(n)+epsZ_temp;
		[ errFlag, vecRho_plus ] = funchRho( rhoArgs, vecZ_temp );
		if (errFlag)
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchRho() failed for +%0.3g to element %d.", epsZ_temp, n ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealvector(vecRho_plus,sizeRho) );
		end
		%
		vecZ_temp = vecZ;
		vecZ_temp(n) = vecZ(n)-epsZ_temp;
		[ errFlag, vecRho_minus ] = funchRho( rhoArgs, vecZ_temp );
		if (errFlag)
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchRho() failed for -%0.3g to element %d.", epsZ_temp, n ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealvector(vecRho_minus,sizeRho) );
		end
		%
		matRhoP(:,n) = ( vecRho_plus - vecRho_minus ) / (2.0*epsZ_temp);
		%
		clear vecRho_plus;
		clear vecRho_minus;
		clear epsZ_temp;
	end
	%
	if ( mygetfield(prm,"useCustomOmega",false) )
		funchOmega   = mygetfield(prm,"funchOmega");
		funchOmegaDRho  = mygetfield(prm,"funchOmegaDRho");
		funchOmegaDRho2 = mygetfield(prm,"funchOmegaDRho2");
	else
		funchOmega   = @(rho)( 0.5 * sum(rho.^2) );
		funchOmegaDRho  = @(rho)( rho );
		funchOmegaDRho2 = @(rho)( eye(sizeRho,sizeRho) );
	end
	omega0      = funchOmega(vecRho0);
	vecOmegaP0  = funchOmegaDRho(vecRho0);
	matOmegaPP0 = funchOmegaDRho2(vecRho0);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega0) );
		assert( omega0 >= 0.0 );
		assert( isrealvector(vecOmegaP0,sizeRho) );
		assert( isrealarray(matOmegaPP0,[sizeRho,sizeRho]) );
	end
	if ( 0.0 == omega0 )
		msg_notify( verbLev, thisFile, __LINE__, "Initial omega is already zero." );
		vecDelta = zeros(size(vecZ));
		retCode = RETCODE__SUCCESS;
		return;
	end
	%
	% For this next part, we assume an orientation for "vec" and "mat"...
	if ( issize(vecRho0,[1,sizeRho]) )
		vecRho0 = vecRho0';
	end
	if ( issize(vecOmegaP0,[1,sizeRho]) )
		vecOmegaP0 = vecOmegaP0';
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecRho0,[sizeRho,1]) );
		assert( isrealarray(vecOmegaP0,[sizeRho,1]) );
	end
	%
	vecG = matRhoP'*vecOmegaP0;
	matH = matRhoP'*matOmegaPP0*matRhoP;
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecG,[sizeZ,1]) );
		assert( isrealarray(matH,[sizeZ,sizeZ]) );
	end
	hScale = sqrt(sum(sum(matH.^2)/sizeZ));
	if ( valLev >= VALLEV__MEDIUM )
		assert( hScale > 0.0 );
		assert( sum(sum( (matH'-matH).^2 )/sizeZ) < eps150*(hScale^2) );
	end
	%
	if ( 0.0 == norm(vecG) )
		msg_notify( verbLev, thisFile, __LINE__, "Initial gradient is zero." );
		% Note: Hessian might not be pos-def,
		%  in which case this might not be a local minimum.
		%  But, such analysis is not supported here.
		vecDelta = zeros(size(vecZ));
		retCode = RETCODE__SUCCESS;
		return;
	end
	%
	%
	%
	% Perform "preliminary regularization" of matH.
	matR = mygetfield( prm, "matR", eps025*hScale*eye(sizeZ,sizeZ) );
	mu0 = mygetfield( prm, "mu0", 1.0E-8 );
	mu1 = mygetfield( prm, "mu1", 1.0E8 );
	muStep = mygetfield( prm, "muStep", 10.0 );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(matR,[sizeZ,sizeZ]) );
		[ foo, errFlag ] = chol(matR);
		if (errFlag)
			msg_error( verbLev, thisFile, __LINE__, ...
			  "Preliminary regularization matrix itself is not positive definite." );
			prelimReguMatIsPosDef = false;
			assert(prelimReguMatIsPosDef);
		end
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
		matH_trial = matH + mu_trial*matR;
		[ foo, errFlag ] = chol(matH_trial);
	end
	matH = matH_trial;
end

	%
	%
	% For BT proper, use Lev or LevMarq.
	%%%msg( thisFile, __LINE__, "TODO: Properly handle backtracking and bounds." );
	vecDelta = -matH\vecG;
	if ( issize(vecZ,[1,sizeZ]) )
		vecDelta = vecDelta';
	end
	
	
	echo__vecG = vecG
	echo__matH = matH
	echo__vecDelta = vecDelta
	
	%
	retCode = RETCODE__SUCCESS;
	return;
end
