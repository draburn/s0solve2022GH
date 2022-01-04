function [ vecX, datOut ] = findMin_bfgs( vecX0, funchOmega, funchG, prm=[] )
	thisFile = "findMin_bfgs";
	%msg( thisFile, __LINE__, "PLACEHOLDER HA~ACK!" );
	%error( "Not implemented." );
	% Later: allow funchOmega providing g and funchG providing omega?
	doChecks = true;
	%
	sizeX = size(vecX0,1);
	omega0 = funchOmega( vecX0 );
	vecG0 = funchG( vecX0 );
	if (doChecks)
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealscalar(omega0) );
		assert( 0.0 <= omega0 );
		assert( isrealarray(vecG0,[sizeX,1]) );
	end
	x0Norm = norm(vecX0);
	g0Norm = norm(vecG0);
	%
	%
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	omegaTol = mygetfield( prm, "omegaTol", eps^0.75 * omega0 + eps  );
	gNormTol = mygetfield( prm, "gNormTol", eps^0.75 * g0Norm + eps * sizeX );
	deltaXNormTol = mygetfield( prm, "deltaXNormTol", eps^0.50 * x0Norm + eps^0.75 * sizeX );
	deltaOmegaTol = mygetfield( prm, "deltaOmegaTol", eps^0.75 * omega0 + eps );
	deltaOmegaRelTol = mygetfield( prm, "deltaOmegaRelTol", eps^0.50 );
	btLimit = mygetfield( prm, "btLimit", 20 );
	if (doChecks)
		assert( isrealscalar(iterLimit) )
		assert( 0 <= iterLimit )
		%
		assert( isrealscalar(omegaTol) );
		assert( isrealscalar(gNormTol) );
		assert( isrealscalar(deltaOmegaTol) );
		assert( isrealscalar(deltaXNormTol) );
		assert( isrealscalar(deltaOmegaRelTol) );
		assert( 0.0 <= omegaTol );
		assert( 0.0 <= gNormTol );
		assert( 0.0 <= deltaOmegaTol );
		assert( 0.0 <= deltaXNormTol );
		assert( 0.0 <= deltaOmegaRelTol );
		%
		assert( isrealscalar(btLimit) )
		assert( 0 <= btLimit )
	end
	%
	%
	vecX = vecX0;
	datOut = [];
	omega = omega0;
	vecG = vecG0;
	iterCount = 0;
	while (1)
		if ( omega <= omegaTol )
			msg( thisFile, __LINE__, "Reached omegaTol." );
			return;
		end
		if ( norm(vecG) <= gNormTol )
			msg( thisFile, __LINE__, "Reached gNormTol." );
			return;
		end
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			return;
		end
		%
		%
		%
		% SELECT STEP DIRECTION vvvvv
		if ( norm(vecG) > 0.01 )
			vecDeltaX = -0.01*vecG/norm(vecG);
		else
			vecDeltaX = -vecG;
		end
		% SELECT STEP DIRECTION ^^^^^
		assert( isrealarray(vecDeltaX,[sizeX,1]) );
		%
		%
		%
		%  DO BACKTRACKING and stuff vvvvv
		vecX_next = vecX + vecDeltaX;
		omega_next = funchOmega( vecX_next );
		if (doChecks)
			assert( isrealscalar(omega_next) );
			assert( 0.0 <= omega_next );
		end
		btCount = 0;
		while ( omega_next >= omega )
			btCount++;
			if ( btCount > btLimit )
				msg( thisFile, __LINE__, "Reached btLimit." );
				return;
			end
			%
			vecDeltaX /= 10.0;
			vecX_next = vecX + vecDeltaX;
			omega_next = funchOmega( vecX_next );
			if (doChecks)
				assert( isrealscalar(omega_next) );
				assert( 0.0 <= omega_next );
			end
		end
		%  DO BACKTRACKING and stuff ^^^^^
		assert( omega_next < omega );
		%
		vecG_next = funchG( vecX_next );
		if (doChecks)
			assert( isrealscalar(omega_next) );
			assert( 0.0 <= omega_next );
			assert( isrealarray(vecG_next,[sizeX,1]) );
		end
		deltaOmega = omega_next - omega;
		%
		%
		% UPDATE MODEL HESSIAN vvvvv
		% UPDATE MODEL HESSIAN ^^^^^
		%
		% Complete iteration.
		vecX = vecX_next;
		omega = omega_next;
		vecG = vecG_next;
		%
		%
		%
		if ( norm(vecDeltaX) <= deltaXNormTol )
			msg( thisFile, __LINE__, "Reached deltaXNormTol." );
			return;
		end
		if ( abs(deltaOmega) <= deltaOmegaTol )
			msg( thisFile, __LINE__, "Reached deltaOmegaTol." );
			return;
		end
		if ( abs(deltaOmega) <= omega*deltaOmegaRelTol )
			msg( thisFile, __LINE__, "Reached deltaOmegaRelTol." );
			return;
		end
	end
return;
end

%!test
%!	clear;
%!	thisFile = "findMin_bfgs test 1";
%!	commondefs;
%!	setprngstates(0);
%!	numFigs = 0;
%!	msg( thisFile, __LINE__, "Yo." );
%!	%
%!	sizeX = 5;
%!	vecXRoot = randn(sizeX,1)
%!	funchOmega = @(x)( 0.5*(norm(x-vecXRoot)^2) );
%!	funchG = @(x)( (x-vecXRoot) );
%!	vecX0 = randn(sizeX,1)
%!	%
%!	tic()
%!	vecXF = findMin_bfgs( vecX0, funchOmega, funchG )
%!	assert( norm(vecXF-vecXRoot) < 1e-2 )
%!	toc()
