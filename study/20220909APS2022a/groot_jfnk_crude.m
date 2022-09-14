function [ vecXBest, fevalCount, matCnvg, datOut ] = groot_jfnk_crude( funchF, vecX0, prm=[] )
	if ( 0 == nargin )
		vecXBest = __FILE__;
		return;
	endif
	%
	groot__commonInit;
	%
	vecXBest = [];
	fevalCount = 0;
	matCnvg = [];
	datOut = [];
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
		if ( fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			break;
		endif
		if ( f <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCCESS: Reached fTol." );
			break;
		endif
		%
		%
		% We hypothetically could include fevalLimit - fevalCount - 1 as a limit inside linsolf.
		funchMatJProd = @(v)(  ( funchF( vecX + (prm.epsFD*v) ) - vecF ) / prm.epsFD  );
		linsolfPrm = mygetfield( prm, "linsolfPrm", [] );
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
		%
		prm.btCoeff = mygetfield( prm, "btCoeff", 0.5 );
		while ( f >= fPrev )
			if ( fevalCount >= prm.fevalLimit )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
				break;
			endif
			if ( norm(vecDelta) <= prm.stepTol )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached stepTol." );
				break;
			endif
			vecDelta *= prm.btCoeff;
			vecX = vecXPrev + vecDelta;
			vecF = funchF( vecX );
			fevalCount++;
			f = norm(vecF);
		endwhile
		%
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
			break;
		elseif ( fPrev - f <= prm.fallTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fallTol on fPrev - f." );
			break;
		endif
	endwhile
	datOut.elapsedTime = time() - startTime;
return;
endfunction
