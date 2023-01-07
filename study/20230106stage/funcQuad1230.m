function [ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW, noisePrm=zeros(3,2) );
	sizeX = size(matX,1);
	sizeS = size(matAS,1);
	sizeW = size(matAW,1);
	numPts = size(matX,2);
	%
	assert( ~isempty(matAS) || ~isempty(matAW) )
	%
	xVarC = noisePrm(1,1);
	fVarC = noisePrm(1,2);
	xVarS = noisePrm(2,1);
	fVarS = noisePrm(2,2);
	xVarW = noisePrm(3,1);
	fVarW = noisePrm(3,2);
	%
	matRXC = xVarC*randn(sizeX,numPts);
	matD = matX - vecXCrit - matRXC;
	if (~isempty(matAS))
		matRXS = 1.0 + xVarS*randn(sizeX,numPts);
		matRFS = 1.0 + fVarS*randn(sizeS,numPts);
		matG = matRXS.*(  matAS' * ((matRFS.^2).*( matAS * (matRXS.*matD) )) );
	endif
	if (~isempty(matAW))
		matRXW = 1.0 + xVarW*randn(sizeX,numPts);
		matRFW = 1.0 + fVarW*randn(sizeW,numPts);
		if (~isempty(matAS))
			matG += matRXW.*(  matAW' * ((matRFW.^2).*( matAW * (matRXW.*matD) )) );
		else
			matG = matRXW.*(  matAW' * ((matRFW.^2).*( matAW * (matRXW.*matD) )) );
		endif
	endif
	matRFC = fVarC*randn(1,numPts);
	rvecF = fCrit + (sum( matD.*matG, 1 )/2.0) + abs(matRFC);
return;
endfunction

%!test
%!	setprngstates(0);
%!	sizeX = 10
%!	numPts = 5
%!	sizeFactor = 0.4
%!	densityFactor = 0.4
%!	%
%!	vecXCrit = randn(sizeX,1);
%!	fCrit = 0.0;
%!	matAS = sprandn(sizeX,sizeX,densityFactor) + diag(randn(sizeX,1));
%!	matAW = randn(round(sizeFactor*sizeX),sizeX);
%!	%
%!	matX = [ randn(sizeX,numPts), vecXCrit ]
%!	%
%!	tic()
%!	[ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW )
%!	toc();

%!test
%!	setprngstates(0);
%!	sizeX = 1E5
%!	numPts = 10
%!	sizeFactor = 2.0E-4
%!	densityFactor = 2.0E-4
%!	%
%!	startTime = time();
%!	vecXCrit = randn(sizeX,1);
%!	fCrit = 0.0;
%!	matAS = sprandn(sizeX,sizeX,densityFactor) + diag(randn(sizeX,1));
%!	matAW = randn(round(sizeFactor*sizeX),sizeX);
%!	compare_mem = [ (sizeof(matAS)+sizeof(matAW)), sizeof(fCrit)*sizeX*sizeX, (sizeof(matAS)+sizeof(matAW))/(sizeof(fCrit)*sizeX*sizeX) ]
%!	matX = [ ones(sizeX,numPts), vecXCrit ];
%!	msg( __FILE__, __LINE__, sprintf( "Generating parameters took %fs.", time()-startTime ) );
%!	%
%!	startTime = time();
%!	[ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW );
%!	msg( __FILE__, __LINE__, sprintf( "Calculating %d point(s) took %fs.", size(matX,2), time()-startTime ) );
%!	rvecF
