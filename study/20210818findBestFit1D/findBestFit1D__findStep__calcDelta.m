function [ vecDelta, retCode, datOut ] = findBestFit1D__findStep__calcDelta( vecZ, omega, vecG, matH, mu, prm )
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
	%
	%
	matZBounds = mygetfield( prm, "matZBounds", [] );
	if (~isempty(matZBounds))
		assert( issize(matZBounds,[sizeZ,2]) );
		for n=1:sizeZ
			assert( isrealorinfscalar(matZBounds(n,1)) );
			assert( isrealorinfscalar(matZBounds(n,2)) );
			assert( matZBounds(n,1) <= vecZ(n) );
			assert( vecZ(n) <= matZBounds(n,2) );
		end
		%
		vecElemIsBounded = zeros(sizeZ,1);
		checkOOB = true;
		while (checkOOB)
			checkOOB = false;
			for n=1:sizeZ
			if (~vecElemIsBounded(n))
				if ( vecZ(n) + vecDelta(n) < matZBounds(n,1) )
					checkOOB = true;
					vecDelta(n) = matZBounds(n,1) - vecZ(n);
					vecElemIsBounded(n) = 1;
				elseif ( vecZ(n) + vecDelta(n) > matZBounds(n,2) )
					checkOOB = true;
					vecDelta(n) = matZBounds(n,2) - vecZ(n);
					vecElemIsBounded(n) = 1;
				end
			end
			end
			if (~checkOOB)
				break;
			end
			%
			matV = [];
			for n=1:sizeZ
			if (~vecElemIsBounded(n))
				vecTemp = zeros(sizeZ,1);
				vecTemp(n) = 1;
				matV = [ matV, vecTemp ];
			end
			end
			%
			matAMod = matV' * matA * matV;
			vecGMod = matV' * ( vecG + matH * vecDelta );
			[ matRMod, errFlag ] = chol( matAMod );
			assert(~errFlag);
			vecDeltaMod = -(matRMod \ ( matRMod' \ vecGMod ));
			vecDelta = vecDelta + matV*vecDeltaMod;
		end
	end
	%
	retCode = RETCODE__SUCCESS;
	return;
end
