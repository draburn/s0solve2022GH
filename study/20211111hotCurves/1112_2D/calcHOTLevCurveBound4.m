function [ matY, vecS ] = calcHOTLevCurveBound4( modelFuncPrm, vecX0, prm=[] )
	thisFile = "calcHOTLevCurveBound4";
	%msg( thisFile, __LINE__, "SIMPLY LIMITING STEP SIZE WON'T WORK; STILL PURELY GRAD.");
	%error( "SIMPLY LIMITING STEP SIZE WON'T WORK; STILL PURELY GRAD." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	maxStepSize = 0.01;
	iterMax = 1000;
	%
	vecX = vecX0;
	matY(:,1) = vecX;
	%
	nextPtIndex = 2;
	nextPtDist = maxStepSize;
	for n=1:iterMax
		%
		%echo__vecX = vecX
		% Get our candiate delta.
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		gNorm = norm(vecG);
		if ( 0.0 == gNorm )
			return;
		end
		vecDelta = -maxStepSize*vecG/norm(vecG);
		%
		newPtDist = norm( vecX + vecDelta - vecX0 );
		if ( 2 == nextPtIndex )
			newDangle = 0.0;
		else
			newDAngle = 1.0 - (vecDelta'*(vecX-vecX0)/( norm(vecDelta)*norm(vecX-vecX0) ));
		end
		%
		if ( newPtDist < 0.99*nextPtDist )
			vecX += vecDelta;
		elseif ( newDangle < 0.1 )
			% We'll expand maxStepSize.
		
		
		if ( newDAngle < 0.1 )
			if ( newPtDist < 0.99*maxStepSize )
			else
			else
			end
		elseif ( newPtDist > maxStepSize )
		else
		end
		%
		matY(:,nextPtIndex) = vecX;
		nextPtIndex++;
		nextPtDist += maxStepSize;
	end
	%
return;
end
