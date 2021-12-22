function vecG = testFunc_gForLSODE( vecX, testFuncPrm )
	assert( isrealarray(vecX,[testFuncPrm.sizeX,1]) );
	%
	vecY = vecX-testFuncPrm.vecXE;
	vecF = testFuncPrm.vecFE + testFuncPrm.matJ*vecY;
	for n=1:testFuncPrm.sizeF
		vecF(n) += 0.5*vecY'*(testFuncPrm.ary3K(:,:,n))*vecY;
	end
	%
	matJ = testFuncPrm.matJ;
	for n=1:testFuncPrm.sizeF
		matJ(n,:) += ( testFuncPrm.ary3K(:,:,n)*vecY )';
	end
	%
	vecG = matJ'*vecF;
	%
return
