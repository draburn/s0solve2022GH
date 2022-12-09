function [ rvecF, matG ] = funcSimpleQuad( matX, vecX0, f0, matH, noiseDat = zeros(2,3) )
	matD = matX - vecX0;
	matD = matD .* ( 1.0 + noiseDat(2,1)*randn(size(matD)) ) + noiseDat(1,1)*randn(size(matD));
	%
	rvecF = f0 + 0.5*sum( matD .* ( matH * matD ), 1 );
	matG = matH * matD;
	%
	rvecF = rvecF .* ( 1.0 + noiseDat(2,2)*randn(size(rvecF)) ) + noiseDat(1,2)*randn(size(rvecF));
	matG = matG .* ( 1.0 + noiseDat(2,3)*randn(size(matG)) ) + noiseDat(1,3)*randn(size(matG));
return;
endfunction

%!test
%!	setprngstates();
%!	sizeX = 5
%!	numPts = 10
%!	%
%!	vecX0 = randn(sizeX,1);
%!	f0 = abs(randn());
%!	matA = randn(sizeX,sizeX);
%!	matH = matA' * matA;
%!	noiseDat = [ 0.0, 0.0, 0.0; 0.0, 0.0, 0.0 ];
%!	%
%!	matX = [ randn(sizeX,numPts), vecX0 ]
%!	%
%!	[ rvecF, matG ] = funcSimpleQuad( matX, vecX0, f0, matH, noiseDat )
