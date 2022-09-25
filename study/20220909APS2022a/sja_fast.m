function [ matJ, datOut ] = sja_fast( matV, matW, prm=[] )
	% AKA sja_omp
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNumNZEPerRow = mygetfield( prm, "maxNumNZEPerRow", sizeK-1 );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	%abortOnBadRow = mygetfield( prm, "abortOnBadRow", true );
	%useVWeight = mygetfield( prm, "useVWeight", false );
	%resumeOnRow = mygetfield( prm, "resumeOnRow", 0 );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeX) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumNZEPerRow) );
		assert( maxNumNZEPerRow <= sizeK );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
		%assert( isscalar(abortOnBadRow) );
		%assert( islogical(abortOnBadRow) );
		%assert( isscalar(useVWeight) );
		%assert( islogical(useVWeight) );
		%assert( 0 == resumeOnRow );
	endif
	%
	matA = matV';
	matE = 1+sparse(sizeF,sizeX); % 1 = not NZE; 0 = is NZE.
	vecProcess = 1+zeros(sizeF,1);
	matNZEList = zeros(sizeF,0);
	matD = matW;
	doMainLoop = true;
	matJ = zeros(sizeF,sizeX);
	for numNZEPerRow = 1 : maxNumNZEPerRow
		vecProcess = sumsq( matD, 2 ) > tol*sumsq( matW, 2 );
		if ( 0 == sum(vecProcess) )
			break;
		endif
		%%%matDAXE = abs(matD * matA) .* matE;
		%%%matE
		[ foo, vecNextNZE ] = max( abs(matD * matA) .* matE, [], 2 );
		matNZEList = [ matNZEList, vecNextNZE ];
		%%%[ sumsq(matD,2), sumsq(matW,2), vecProcess ](1,:)
		%%%matNZEList(1,:)
		for nf = (1 : sizeF)(vecProcess)
			matE(nf,matNZEList(nf,end)) = 0;
			%%%matJ(nf,unique(matNZEList(nf,:))) = matA(:,unique(matNZEList(nf,:))) \ (matW(nf,:)');
			matJ(nf,matNZEList(nf,:)) = matA(:,matNZEList(nf,:)) \ (matW(nf,:)');
		endfor
		matD = matW - matJ*matV;
	endfor
	%
	datOut = [];
	return;
return;
endfunction
