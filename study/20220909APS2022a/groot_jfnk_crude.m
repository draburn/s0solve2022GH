function [ vecXBest, fevalCount, matCnvg, datOut ] = groot_jfnk_crude( funchF, vecX0, fTol=1.0e-6, fallTol=1.0e-7, fevalLimit=1000, prm=[] )
	mydefs;
	vecXBest = [];
	fevalCount = 0;
	matCnvg = [];
	datOut = [];
	prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__FLAGGED );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(fTol) );
	assert( 0.0 < fTol );
	assert( isrealscalar(fallTol) );
	assert( 0.0 < fallTol );
	assert( isposintscalar(fevalLimit) );
	%
	epsFD = mygetfield( prm, "epsFD", eps^(1.0/3.0) );
	assert( isrealscalar(epsFD) );
	assert( 0.0 ~= epsFD );
	%
	vecF0 = funchF( vecX0 );
	fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	f0 = norm(vecF0);
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	f = f0;
	iterCount = 0;
	matCnvg(iterCount+1,:) = [ fevalCount, norm(vecF0) ];
	fBest = f0;
	vecXBest = vecX0;
	while (1)
		if ( fevalCount >= fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			return;
		endif
		if ( f <= fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCCESS: Reached fTol." );
			return;
		endif
		%
		%
		% We hypothetically could include fevalLimit - fevalCount - 1 as a limit inside linsolf.
		funchMatJProd = @(v)(  ( funchF( vecX + (epsFD*v) ) - vecF ) / epsFD  );
		linsolfPrm = [];
		[ vecDelta, linsolfDatOut ] = linsolf( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		fevalCount += linsolfDatOut.fevalCount;
		%
		%
		vecXPrev = vecX;
		fPrev = f;
		%
		vecX = vecXPrev + vecDelta;
		vecF = funchF( vecX );
		fevalCount++;
		f = norm(vecF);
		%
		iterCount++;
		matCnvg(iterCount+1,:) = [ fevalCount, norm(vecF) ];
		if ( f < fBest )
			fBest = f;
			vecXBest = vecX;
		endif
		%
		%
		if ( f > fPrev )
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: f > fPrev." );
			return;
		elseif ( fPrev - f <= fallTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fallTol on fPrev - f." );
			return;
		endif
	endwhile
return;
endfunction
