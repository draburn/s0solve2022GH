function [ rvecF, matG ] = funcQuadCustom1229( matX, vecXCrit, fCrit=0.0, matASparse=[], matAFull=[], prm=[] );
	sz = size(matX,1);
	np = size(matX,2);
	if ( mygetfield( prm, "debugMode", false ) )
		assert( isposintscalar(sz) );
		assert( isposintscalar(np) );
		assert( isrealarray(matX,[sz,np]) );
		assert( isrealarray(vecXCrit,[sz,1]) );
		assert( isrealscalar(fCrit) );
		if (~isempty(matASparse))
			assert( isrealarray(matASparse) );
			assert( size(matASparse,2) == sz );
		endif
		if (~isempty(matAFull))
			assert( isrealarray(matAFull) );
			assert( size(matAFull,2) == sz );
		endif
	endif
	%
	vecXCritVar = mygetfield( prm, "vecXCritVar", 1.0e-8+(1.0e-6*vecXCrit) );
	fCritVar = mygetfield( prm, "fCritVar", 0.0 );
	matASparseVar = mygetfield( prm, "matASparseVar", [] );
	if (isempty(matASparseVar))
		matASparseVar = 1.0e-2*matASparse;
		matASparseVar(find(matASparse)) += 1.0e-4;
	endif
	matAFullVar = mygetfield( prm, "matAFullVar", [] );
	if (isempty(matAFullVar))
		matAFullVar = 1.0e-4+(1.0e-2*matAFull);
	endif
	bSparseVar = mygetfield( prm, "bSparseVar", 1.0e-4 );
	bDensityMax = mygetfield( prm, "bDensityMax", 1.0e-2 );
	bFullVar = mygetfield( prm, "bFullVar", 1.0e-4 );
	bSizeMax = mygetfield( prm, "bSizeMax", 1.0e-2*sz );
	%
	matXCrit = vecXCrit + vecXCritVar.*randn(sz,np);
	rvecFCrit = fCrit + fCritVar*abs(randn(1,np));
	matD = matX - matXCrit;
	matG = zeros(sz,np);
	if ( ~isempty(matASparse) )
		if ( max(abs(matASparseVar)) > 0.0 )
			indexList = find(matASparse);
			matASparse(indexList) += matASparseVar(indexList) .* randn(size(indexList));
		endif
		matG += matASparse' * ( matASparse * matD );
	endif
	if ( ~isempty(matAFull) )
		if ( max(abs(matAFullVar)) > 0.0 )
			matAFull += matAFullVar .* randn(size(matAFull));
		endif
		matG += matAFull' * ( matAFull * matD );
	endif
	if ( 0.0 < bSparseVar && 0.0 < bDensityMax )
		matBSparse = bSparseVar * sprandn( sz, sz, rand()*bDensityMax );
		matG += matBSparse' * ( matBSparse * matD );
	endif
	if ( 0.0 < bFullVar && 0 < bSizeMax )
		matBFull = bFullVar * randn( floor(rand()*bSizeMax), sz );
		matG += matBFull' * ( matBFull * matD );
	endif
	matG = full( matG );
	rvecF = rvecFCrit + sum( matD .* matG, 1 )/2.0;
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
%!	fCrit = abs(randn());
%!	matASparse = sprandn(sizeX,sizeX,densityFactor) + diag(randn(sizeX,1));
%!	matAFull = randn(round(sizeFactor*sizeX),sizeX);
%!	%
%!	matX = zeros(sizeX,numPts);
%!	%
%!	prm = [];
%!	[ rvecF, matG ] = funcQuadCustom1229( matX, vecXCrit, fCrit, matASparse, matAFull, prm )

%!test
%!	setprngstates(0);
%!	sizeX = 2E4
%!	numPts = 1
%!	sizeFactor = 1.0E-2
%!	densityFactor = 1.0E-3
%!	%
%!	tic();
%!	msgnnl( __FILE__, __LINE__, "Generating data... " );
%!	vecXCrit = randn(sizeX,1);
%!	fCrit = abs(randn());
%!	matASparse = sprandn(sizeX,sizeX,densityFactor); + diag(randn(sizeX,1));
%!	matAFull = randn(round(sizeFactor*sizeX),sizeX);
%!	%
%!	matX = zeros(sizeX,numPts);
%!	%
%!	prm = [];
%!	prm.vecXCritVar = 1.0E-8+(1.0E-6*vecXCrit);
%!	prm.fCritVar = 0.0;
%!	%prm.matASparseVar = 1.0E-2*matASparse;
%!	%prm.matASparseVar(find(matASparse)) += 1.0E-4;
%!	prm.matASparseVar = 0.0;
%!	%prm.matAFullVar = 1.0E-4+(1.0E-2*matAFull);
%!	prm.matAFullVar = 0.0;
%!	prm.bSparseVar = 0.0;
%!	prm.bDensityMax = 1.0E-7;
%!	prm.bFullVar = 1.0E-2;
%!	prm.bSizeMax = 1.0E-3*sizeX;
%!	toc();
%!	%
%!	tic();
%!	msgnnl( __FILE__, __LINE__, "Evaluating point... " );
%!	[ rvecF, matG ] = funcQuadCustom1229( matX, vecXCrit, fCrit, matASparse, matAFull, prm );
%!	toc();
