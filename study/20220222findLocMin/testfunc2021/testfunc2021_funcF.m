function [ vecFVals, matJVals ] = testfunc2021_funcF( vecXVals, testFuncPrm )
	numPts = size(vecXVals,2);
	assert( isposintscalar(numPts) );
	assert( isrealarray(vecXVals,[testFuncPrm.sizeX,numPts]) );
	%
	if ( 1==numPts)
		vecY = vecXVals-testFuncPrm.vecXE;
		vecFVals = testFuncPrm.vecFE + testFuncPrm.matJ*vecY;
		for n=1:testFuncPrm.sizeF
			vecFVals(n) += 0.5*vecY'*(testFuncPrm.ary3K(:,:,n))*vecY;
		end
		if ( 1 == nargout )
			return;
		end
		matJVals = testFuncPrm.matJ;
		for n=1:testFuncPrm.sizeF
			matJVals(n,:) += ( testFuncPrm.ary3K(:,:,n)*vecY )';
		end
	else
		vecYVals = vecXVals - testFuncPrm.vecXE;
		vecFVals = testFuncPrm.vecFE + testFuncPrm.matJ*vecYVals;
		for n=1:testFuncPrm.sizeF
			vecFVals(n,:) += 0.5*sum( vecYVals.*( testFuncPrm.ary3K(:,:,n)*vecYVals ), 1 );
		end
		if ( 1 == nargout )
			return;
		end
		error( "Calculation of matJVals for numPts > 1 is not supported." );
	end
return
