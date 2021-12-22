function [ vecF, matJ, vecG, matH ] = func1111_evalDeriv( vecX, funcPrm )
	assert( isrealarray(vecX,[funcPrm.sizeX,1]) );
	%
	vecY = vecX-funcPrm.vecXE;
	vecF = funcPrm.vecFE + funcPrm.matJ*vecY;
	for n=1:funcPrm.sizeF
		vecF(n) += 0.5*vecY'*(funcPrm.ary3K(:,:,n))*vecY;
	end
	%
	matJ = funcPrm.matJ;
	for n=1:funcPrm.sizeF
		matJ(n,:) += ( funcPrm.ary3K(:,:,n)*vecY )';
	end
	%
	vecG = matJ'*vecF;
	%
	matH = matJ'*matJ;
	for n=1:funcPrm.sizeF
		matH += vecF(n)*funcPrm.ary3K(:,:,n);
	end
	%
return
