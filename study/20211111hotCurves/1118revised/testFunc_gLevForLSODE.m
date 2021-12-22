function vecGLev = testFunc_gLevForLSODE( vecX, testFuncPrm, vecX0, barrierDistance )
	%echo__testFuncPrm_sizeX = testFuncPrm.sizeX
	%echo__size_vecX_1 = size(vecX,1)
	%echo__size_vecX_2 = size(vecX,2)
	%echo__vecX = vecX
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
	z = norm( vecX - vecX0 );
	if ( z <= barrierDistance )
		vecGLev = vecG;
		return;
	end
	%echo__z = z
	%echo__bd = barrierDistance
	%%%vecGLev = vecG + (vecX-vecX0) * ( exp( (z/barrierDistance)^2 ) - 1.0 );
	%%%vecGLev = vecG + (vecX-vecX0) * ( ((z/barrierDistance)^4) - 1.0 );
	vecGLev = vecG + (vecX-vecX0) * ( exp(( (z/barrierDistance) - 1.0 )^2) - 1.0 );
	%
return
