function [ vecDelta, retCode, datOut ] = findBestFit1D__findStep__calcDelta( omega, vecG, matH, mu, prm )
	%
	% Init
	commondefs;
	thisFile = "findBestFit1D__findStep__calcDelta";
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
	sizeZ = size(vecG,1);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega) );
		assert( 0.0 <= omega );
		assert( isrealarray(vecG,[sizeZ,1]) );
		assert( isrealarray(matH,[sizeZ,sizeZ]) );
	end
	%
	hScale = sqrt(sum(sum(matH.^2)/sizeZ));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(hScale) );
		assert( 0.0 < hScale );
	end
	%
	if ( mygetfield(prm,"useLevMarq",true) )
		matD = diag(diag(matH));
	else
		matD = hScale*eye(sizeZ,sizeZ);
	end
	%
	matA = matH + (mu*matD);
	[ matR, errFlag ] = chol( matA );
	if ( errFlag )
		msg_error( verbLev, thisFile, __LINE__, "ERROR: Cholesky factorization failed." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	%
	vecDelta = -(matR \ ( matR' \ vecG ));
	%
	msg( thisFile, __LINE__, "TODO: Add loop to handle boundaries." );
	if (0)
	vecBounded = zeros(sizeZ,1);
	checkOOB = true;
	while (checkOOB)
		checkOOB = false;
		for n=1:sizeZ
		if (~vecBounded(n))
			if ( vecZ(n) + vecDelta(n) < vecZMin(n) )
				checkOOB = true;
				vecDelta(n) = vecZMin(n) - vecZ(n);
				vecBounded(n) = 1;
			elseif ( vecZ(n) + vecDelta(n) > vecZMax(n) )
				checkOOB = true;
				vecDelta(n) = vecZMax(n) - vecZ(n);
				vecBounded(n) = 1;
			end
		end
		end
		%
		msg( thisFile, __LINE__, "TODO: Properly update vecDelta when bounded." );
		%
	end
	end
	%
	retCode = RETCODE__SUCCESS;
	return;
end
