function [ rvecF, matG ] = funcQuadWNoiseSparse1229( matX, vecX0, f0, matH, noiseDat=zeros(4,2), noiseDensityMax=[] )
	matD = matX - vecX0;
	matD = matD .* ( 1.0 + noiseDat(1,1)*randn(size(matD)) ) + noiseDat(1,2)*randn(size(matD));
	%
	if ( isempty(noiseDensityMax) )
		noiseDensityMax = nnz(matH)/prod(size(matH));
	endif
	matH = matH .* ( 1.0 + noiseDat(2,1)*randn(size(matH)) ) + noiseDat(2,2)*sprandn( size(matH,1), size(matH,2), rand()*noiseDensityMax );
	matH = (matH'+matH)/2.0;
	rvecF = f0 + sum( matD .* ( matH * matD ), 1 )/2.0;
	matG = matH * matD;
	%
	rvecF = rvecF .* ( 1.0 + noiseDat(3,1)*randn(size(rvecF)) ) + noiseDat(3,2)*randn(size(rvecF));
	matG = matG .* ( 1.0 + noiseDat(4,1)*randn(size(matG)) ) + noiseDat(4,2)*randn(size(matG));
return;
endfunction

%!test
%!	setprngstates();
%!	sizeX = 5
%!	numPts = 10
%!	densityFactor = 0.4
%!	%
%!	vecX0 = randn(sizeX,1);
%!	f0 = abs(randn());
%!	matA = sprandn( sizeX, sizeX, densityFactor );
%!	matH = matA' * matA;
%!	noiseDat = sqrt(eps)+zeros(4,2)
%!	%
%!	matX = [ randn(sizeX,numPts), vecX0 ]
%!	%
%!	[ rvecF, matG ] = funcQuadWNoiseSparse1229( matX, vecX0, f0, matH, noiseDat )
