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
	% Note that this may lead to the same "corner" vecDelta being used
	% for multiple values of mu.
	% A straightforward optimization would be to avoid that.
	matZBounds = mygetfield( prm, "matZBounds", [] );
	if (~isempty(matZBounds))
		assert( issize(matZBounds,[sizeZ,2]) );
		for n=1:sizeZ
			assert( isrealorinfscalar(matZBounds(n,1)) );
			assert( isrealorinfscalar(matZBounds(n,2)) );
			%assert( matZBounds(n,1) <= vecZ(n) );
			%assert( vecZ(n) <= matZBounds(n,2) );
			assert( matZBounds(n,1)-eps075*abs(matZBounds(n,1)) <= vecZ(n)+eps075*abs(vecZ(n)) );
			assert( vecZ(n)-eps075*abs(vecZ(n)) <= matZBounds(n,2)+eps075*abs(matZBounds(n,2)) );
		end
		%
		vecElemIsBounded = zeros(sizeZ,1);
		while (true)
			% Check boundaries; identify closes one per straight-line.
			nOfClosestBound = 0; % 0 indicates no bound.
			cOfClosestBound = 2.0; % Largest meaningful value is 1.0.
			for n=1:sizeZ
			if (~vecElemIsBounded(n))
				c = 3.0; % Largest meaningful value is 1.0 (or 2.0 until set).
				if ( 0.0 == vecDelta(n) )
					if (  (vecZ(n) < matZBounds(n,1))  ||  (vecZ(n) > matZBounds(n,2)) )
						c = 0.0;
						% Not sure if this ever comes up,
						% but, we're going to divide by vecDelta(n) below.
					end
				elseif ( vecZ(n) + vecDelta(n) < matZBounds(n,1) )
					c = (matZBounds(n,1) - vecZ(n)) / vecDelta(n);
				elseif ( vecZ(n) + vecDelta(n) > matZBounds(n,2) )
					c = (matZBounds(n,2) - vecZ(n)) / vecDelta(n);
				end
				if ( c < cOfClosestBound )
					cOfClosestBound = c;
					nOfClosestBound = n;
				end
			end
			end
			if ( valLev >= VALLEV__MEDIUM )
				assert( 0.0 <= cOfClosestBound );
				assert( 0 <= nOfClosestBound );
				assert( nOfClosestBound <= sizeZ );
			end
			%
			if ( 0 == nOfClosestBound )
				% Haven't crossed any boundaries.
				break;
			end
			if ( valLev >= VALLEV__MEDIUM )
				assert( cOfClosestBound <= 1.0 );
			end
			%
			% Apply bound.
			n = nOfClosestBound;
			if ( vecZ(n) + vecDelta(n) < matZBounds(n,1) )
				vecDelta(n) = matZBounds(n,1) - vecZ(n);
			elseif ( vecZ(n) + vecDelta(n) > matZBounds(n,2) )
				vecDelta(n) = matZBounds(n,2) - vecZ(n);
			else
				msg_error( verbLev, thisFile, __LINE__, "INTERAL ERROR." );
				retCode = RETCODE__INTERNAL_INCONSISTENCY;
				return;
			end
			vecElemIsBounded(n) = 1;
			%
			% Optimize all other elements.
			matV = [];
			for n=1:sizeZ
			if (~vecElemIsBounded(n))
				vecTemp = zeros(sizeZ,1);
				vecTemp(n) = 1;
				matV = [ matV, vecTemp ];
			end
			end
			if (isempty(matV))
				% We've hit a corner, nothing left to optimize.
				break;
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
	%
	%
	retCode = RETCODE__SUCCESS;
	return;
end
