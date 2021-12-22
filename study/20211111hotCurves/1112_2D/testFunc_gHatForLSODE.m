function [ vecGHat ] = testFunc_gHatForLSODE( vecX, testFuncPrm )
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
	fNorm = norm(vecF);
	jNorm = sqrt(sum(sum(matJ'*matJ)))/testFuncPrm.sizeX;
	vecGHat = vecG/( norm(vecG) + (0.1*fNorm*jNorm) );
	% NOT ACTUALLY G HAT.
	% LSODE BONKS WHEN G HAT GOES TO ZERO.
	%
return
