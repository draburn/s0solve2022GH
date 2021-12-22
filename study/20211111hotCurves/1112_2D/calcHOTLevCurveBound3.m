function [ matY, vecS ] = calcHOTLevCurveBound3( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTLevCurveBound3";
	msg( thisFile, __LINE__, "SIMPLY LIMITING STEP SIZE WON'T WORK; STILL PURELY GRAD.");
	error( "SIMPLY LIMITING STEP SIZE WON'T WORK; STILL PURELY GRAD." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	maxStepSize = 0.01;
	iterMax = 10000;
	%
	vecX = vecX0;
	matY(:,1) = vecX;
	vecS(1) = 0.0;
	%
	nextPtIndex = 2;
	nextPtDist = maxStepSize;
	for n=1:iterMax
		%
		% Get our candiate delta.
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		gNorm = norm(vecG);
		if ( 0.0 == gNorm )
			return;
		end
		vecDelta = -maxStepSize*vecG/norm(vecG);
		%
		if ( norm( vecX + vecDelta - vecX0 ) < nextPtDist )
			% Just take the delta.
			vecX += vecDelta;
		else
			% This puts us past nextPtDist.
			% We'll find the intersection and take it.
			currPtDist = norm(vecX-vecX0);
			assert( currPtDist <= nextPtDist );
			deltaNorm = norm(vecDelta);
			crossTerm = vecDelta'*(vecX-vecX0);
			underSqrtTerm = (crossTerm^2) + (deltaNorm^2)*( (nextPtDist^2) - (currPtDist^2) );
			s = ( sqrt(underSqrtTerm) - crossTerm ) / (deltaNorm^2);
			assert( 0.0 < s );
			assert( s <= 1.0 );
			vecDelta *= s;
			vecX += vecDelta;
			%
			matY(:,nextPtIndex) = vecX;
			vecS(nextPtIndex) = norm(vecX-vecX0);
			nextPtIndex++;
			nextPtDist += 0.3*maxStepSize;
		end
	end
	%
return;
end
