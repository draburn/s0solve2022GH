function [ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW, prm=[] );
	sizeX = size(matX,1);
	sizeS = size(matAS,1);
	sizeW = size(matAW,1);
	numPts = size(matX,2);
	%
	assert( ~isempty(matAS) || ~isempty(matAW) )
	%
	prm.xVarC = mygetfield( prm, "xVarC", 0.0 );
	prm.fVarC = mygetfield( prm, "fVarC", 0.0 );
	prm.xVarS = mygetfield( prm, "xVarS", 1.0E-2 );
	prm.fVarS = mygetfield( prm, "fVarS", 1.0E-2 );
	prm.xVarW = mygetfield( prm, "xVarW", 1.0E-2 );
	prm.fVarW = mygetfield( prm, "fVarW", 1.0E-2 );
	%
	matRXC = prm.xVarC*randn(sizeX,numPts);
	matD = matX - vecXCrit - matRXC;
	if (~isempty(matAS))
		matRXS = 1.0 + prm.xVarS*randn(sizeX,numPts);
		matRFS = 1.0 + prm.fVarS*randn(sizeS,numPts);
		matG = matRXS.*(  matAS' * ((matRFS.^2).*( matAS * (matRXS.*matD) )) );
	endif
	if (~isempty(matAW))
		matRXW = 1.0 + prm.xVarW*randn(sizeX,numPts);
		matRFW = 1.0 + prm.fVarW*randn(sizeW,numPts);
		if (~isempty(matAS))
			matG += matRXW.*(  matAW' * ((matRFW.^2).*( matAW * (matRXW.*matD) )) );
		else
			matG = matRXW.*(  matAW' * ((matRFW.^2).*( matAW * (matRXW.*matD) )) );
		endif
	endif
	matRFC = prm.fVarC*randn(1,numPts);
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
%!	prm = [];
%!	tic()
%!	[ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW, prm )
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
%!	prm = [];
%!	startTime = time();
%!	[ rvecF, matG ] = funcQuad1230( matX, vecXCrit, fCrit, matAS, matAW, prm );
%!	msg( __FILE__, __LINE__, sprintf( "Calculating %d point(s) took %fs.", size(matX,2), time()-startTime ) );
%!	rvecF
