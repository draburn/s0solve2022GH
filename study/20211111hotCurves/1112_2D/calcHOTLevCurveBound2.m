function [ matY, vecS ] = calcHOTLevCurveBound2( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTLevCurveBound2";
	msg( thisFile, __LINE__, "DOES NOT WORK.");
	error( "DOES NOT WORK." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	maxStepSize = 0.01;
	iterMax = 10000;
	%
	vecX = vecX0;
	matY(:,1) = vecX;
	vecS(1) = 0.0;
	ptIndex = 2;
	totalStepSize = maxStepSize
	for n=1:iterMax
		echo__vecX = vecX
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		gNorm = norm(vecG);
		if ( 0.0 == gNorm )
			return;
		end
		vecDelta = -maxStepSize*vecG/norm(vecG);
		vecDPre = vecX - vecX0;
		vecDTot = vecDPre + vecDelta;
		if ( norm(vecDTot) >= totalStepSize ) % Have we gone past current boundary?
			dPreNorm = norm(vecDPre)
			deltaNorm = norm(vecDelta);
			crossTerm = vecDPre'*vecDelta;
			inSqrtTerm = (crossTerm^2) + (deltaNorm^2)*((totalStepSize^2) - (dPreNorm^2));
			assert( inSqrtTerm >= 0.0 );
			msg( thisFile, __LINE__, "HERE." );
			if ( inSqrtTerm < (crossTerm^2) ) % Are we (after binding) close to boundary?
				msg( thisFile, __LINE__, "HERE." );
				matY(:,ptIndex) = vecX;
				vecS(ptIndex) = norm(vecDPre);
				ptIndex++;
				totalStepSize += maxStepSize
				%
				% But, make sure we haven't gone past NEXT boundary!
				% So, don't assign next vecX quite yet.
				dPreNorm = norm(vecDPre)
				deltaNorm = norm(vecDelta);
				crossTerm = vecDPre'*vecDelta;
				inSqrtTerm = (crossTerm^2) + (deltaNorm^2)*((totalStepSize^2) - (dPreNorm^2));
				assert( inSqrtTerm >= 0.0 );
			end
			%
			s = ( sqrt(inSqrtTerm) - crossTerm ) / (deltaNorm^2)
			assert( 0.0 < s );
			assert( s <= 1.0 );
			vecDelta *= s;
			vecX += vecDelta;
		else
			vecX += vecDelta;
			msg( thisFile, __LINE__, "You go, boy!" );
		end
	end
	%
return;
end
