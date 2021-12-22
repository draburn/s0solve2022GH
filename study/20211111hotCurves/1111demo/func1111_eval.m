function matF = func1111_eval( matX, funcPrm )
	numPts = size(matX,2);
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[funcPrm.sizeX,numPts]) );
	%
	if ( 1==numPts)
		vecX = matX;
		vecY = vecX-funcPrm.vecXE;
		vecF = funcPrm.vecFE + funcPrm.matJ*vecY;
		for n=1:funcPrm.sizeF
			vecF(n) += 0.5*vecY'*(funcPrm.ary3K(:,:,n))*vecY;
		end
		matF = vecF;
		return;
	else
		matY = matX - repmat( funcPrm.vecXE, [1,numPts] );
		matF = repmat( funcPrm.vecFE, [1,numPts] ) + funcPrm.matJ*matY;
		for n=1:funcPrm.sizeF
			matF(n,:) += 0.5*sum( matY.*( funcPrm.ary3K(:,:,n)*matY ), 1 );
		end
		return
	end
return
