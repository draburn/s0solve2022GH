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
		assert( isrealvector(vecZ,sizeZ) );
	end
	%
	%
	%
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
		funchOmegaP  = mygetfield(prm,"funchOmegaP");
		funchOmegaPP = mygetfield(prm,"funchOmegaPP");
	else
		funchOmega   = @(rho)( 0.5 * sum(rho.^2) );
		funchOmegaP  = @(rho)( rho );
		funchOmegaPP = @(rho)( eye(sizeRho,sizeRho) );
	end
	omega0      = funchOmega(vecRho0)
	vecOmegaP0  = funchOmegaP(vecRho0)
	matOmegaPP0 = funchOmegaPP(vecRho0)
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega0) );
		assert( omega0 >= 0.0 );
		assert( isrealvector(vecOmegaP0,sizeRho) );
		assert( isrealarray(matOmegaPP0,[sizeRho,sizeRho]) );
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
	%
	%
	%
	msg( thisFile, __LINE__, "TODO: Properly handle regularization, backtracking, and bounds." );
	vecDelta = -matH\vecG;
	if ( issize(vecZ,[1,sizeZ]) )
		vecDelta = vecDelta';
	end
	%
	retCode = RETCODE__SUCCESS;
	return;
end
