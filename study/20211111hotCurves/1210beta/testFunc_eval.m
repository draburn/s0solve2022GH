function matF = testFunc_eval( matX, testFuncPrm )
	numPts = size(matX,2);
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[testFuncPrm.sizeX,numPts]) );
	%
	if ( 1==numPts)
		vecX = matX;
		vecY = vecX-testFuncPrm.vecXE;
		vecF = testFuncPrm.vecFE + testFuncPrm.matJ*vecY;
		for n=1:testFuncPrm.sizeF
			vecF(n) += 0.5*vecY'*(testFuncPrm.ary3K(:,:,n))*vecY;
		end
		matF = vecF;
		return;
	else
		matY = matX - repmat( testFuncPrm.vecXE, [1,numPts] );
		matF = repmat( testFuncPrm.vecFE, [1,numPts] ) + testFuncPrm.matJ*matY;
		for n=1:testFuncPrm.sizeF
			matF(n,:) += 0.5*sum( matY.*( testFuncPrm.ary3K(:,:,n)*matY ), 1 );
		end
		return
	end
return
