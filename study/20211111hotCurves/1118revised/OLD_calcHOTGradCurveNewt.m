function matY = OLD_calcHOTGradCurveNewt( modelFuncPrm, vecX0, prm=[] )
	thisFile = "OLD_calcHOTGradCurveNewt";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	maxNumIter = mygetfield( prm, "maxNumIter", 1000 );
	assert( isposintscalar(maxNumIter) );
	%
	vecX = vecX0;
	n = 1;
	while (1)
		matY(:,n) = vecX;
		n++;
		if ( n>maxNumIter )
			return;
		end
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		vecX -= matH\vecG;
		%vecX -= (matJ'*matJ)\(matJ'*vecF);
		if ( norm(vecG) == 0.0 )
			break;
		end
	end
	return;
end
