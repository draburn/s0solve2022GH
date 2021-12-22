function matX = calcGradCurve_chunky( funchG, vecX0, prm=[] )
	thisFile = "calcGradCurve_chunky";
	targetStepSize = 0.2;
	stepsPerChunk = 20;
	gradNormTol = 1e-6;
	iterLimit = 100;
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0) );
	%%%funchXDot_lsode = @(x,t)( -funchG(x) );
	%funchGHatIsh = @(g)( g/(0.1*targetStepSize + norm(g)) );
	funchGHatIsh = @(g)( g/sqrt( targetStepSize^2 + g'*g ) );
	funchXDot_lsode = @(x,t)( -funchGHatIsh(funchG(x)) );
	%
	numPts = 1;
	matX(:,1) = vecX0;
	vecX = vecX0;
	iterCount = 0;
	while (1)
		%if ( omega <= omegaTol )
		%	msg( thisFile, __LINE__, "Reached omegaTol." );
		%	return;
		%end
		vecG = funchG(vecX);
		gNorm = norm(vecG);
		if ( gNorm <= gradNormTol )
			msg( thisFile, __LINE__, "Reached gradNormTol." );
			return;
		end
		if ( iterCount >= iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			return;
		end
		iterCount++;
		%
		deltaT = stepsPerChunk*targetStepSize/norm(funchXDot_lsode(vecX));
		vecT = deltaT*linspace(0.0,1.0,stepsPerChunk+1);
		assert( isrealarray(vecT) );
		matX_new = lsode( funchXDot_lsode, vecX, vecT )';
		assert( isrealarray(matX_new,[sizeX,stepsPerChunk+1]) );
		%
		% Trim any redundant points.
		% This tends to be necessary on the last chunk.
		numPts_new = stepsPerChunk;
		while (1)
			vec1 = matX_new(:,numPts_new);
			vec2 = matX_new(:,numPts_new+1);
			if ( norm( vec1 - vec2 ) < (eps^0.75)*( norm(vec1) + norm(vec2) ) )
				%msg( thisFile, __LINE__, sprintf( "Dropping point %d.", numPts_new+1 ) );
				numPts_new--;
			else
				break;
			end
			if ( numPts_new <= 0 )
				msg( thisFile, __LINE__, "WARNING: NO new points were added." );
				return;
			end
		end
		%
		matX(:,numPts+1:numPts+numPts_new) = matX_new(:,2:numPts_new+1);
		vecX = matX_new(:,end);
		numPts += numPts_new;
	end
	%
return;
end
