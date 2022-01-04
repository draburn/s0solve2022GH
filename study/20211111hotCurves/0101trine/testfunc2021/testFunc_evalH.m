function matH = testFunc_evalH( vecX, testFuncPrm )
	assert( isrealarray(vecX,[testFuncPrm.sizeX,1]) );
	%
	vecY = vecX-testFuncPrm.vecXE;
	vecF = testFuncPrm.vecFE + testFuncPrm.matJ*vecY;
	for n=1:testFuncPrm.sizeF
		vecF(n) += 0.5*vecY'*(testFuncPrm.ary3K(:,:,n))*vecY;
	end
	%
	%omega = (vecF'*vecF)/2.0;
	%
	matJ = testFuncPrm.matJ;
	for n=1:testFuncPrm.sizeF
		matJ(n,:) += ( testFuncPrm.ary3K(:,:,n)*vecY )';
	end
	%
	%vecG = matJ'*vecF;
	%
	matH = matJ'*matJ;
	for n=1:testFuncPrm.sizeF
		matH += vecF(n)*testFuncPrm.ary3K(:,:,n);
	end
	%
return
