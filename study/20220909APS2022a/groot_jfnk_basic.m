function [ vecXBest, grootFlag, fevalCount, datOut ] = groot_jfnk_basic( funchF, vecX0, prm=[] )
	if ( 0 == nargin )
		vecXBest = __FILE__;
		return;
	elseif ( nargin < 2 )
		error( "Too few input arguments." );
	elseif ( 3 < nargin )
		error( "Too many input arguments." );
	elseif ( 4 < nargout )
		error( "Too many output arguments." );
	endif
	groot__commonInit;
	vecXBest = [];
	grootFlag = GROOT_FLAG__VOID;
	fevalCount = 0;
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
	matRecordX(:,iterCount+1) = vecX;
	matInfoA(iterCount+1,:) = [ iterCount, fevalCount, norm(vecX), norm(vecF) ];
	matInfoB = [];
	fBest = f0;
	vecXBest = vecX0;
	doMainLoop = true;
	while (doMainLoop)
		if ( f <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCCESS: Reached fTol." );
			grootFlag = GROOT_FLAG__CNVG;
			doMainLoop = false; % Redundant.
			break;
		elseif ( fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		elseif ( stopsignalpresent() )
			msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		endif
		%
		%
		% We hypothetically could include fevalLimit - fevalCount - 1 as a limit inside linsolf.
		funchMatJProd = @(v)(  ( funchF( vecX + (prm.epsFD*v) ) - vecF ) / prm.epsFD  );
		linsolfPrm = mygetfield( prm, "linsolfPrm", [] );
		[ vecDelta0, linsolfDatOut ] = linsolf( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		fevalCount += linsolfDatOut.fevalCount;
		if ( fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		endif
		%
		%
		vecXPrev = vecX;
		vecFPrev = vecF;
		fPrev = f;
		%
		vecDelta = vecDelta0;
		vecX = vecXPrev + vecDelta;
		vecF = funchF( vecX );
		fevalCount++;
		f = norm(vecF);
		%
		%
		needBT = ( f >= fPrev - prm.fallTol );
		prm.btCoeff = mygetfield( prm, "btCoeff", 0.5 );
		if ( needBT && prm.btCoeff <= 0.0 );
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Reached fallTol with no BT." );
			grootFlag = GROOT_FLAG__FAIL;
			doMainLoop = false; % Redundant.
			break;
		endif
		while ( needBT )
			if ( norm(vecDelta) <= prm.stepTol )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Reached stepTol." );
				grootFlag = GROOT_FLAG__FAIL;
				doMainLoop = false;
				break;
			elseif ( fevalCount >= prm.fevalLimit )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
				grootFlag = GROOT_FLAG__STOP;
				doMainLoop = false;
				break;
			elseif ( stopsignalpresent() )
				msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
				grootFlag = GROOT_FLAG__STOP;
				doMainLoop = false;
				break;
			endif
			vecDelta *= prm.btCoeff;
			vecX = vecXPrev + vecDelta;
			vecF = funchF( vecX );
			fevalCount++;
			f = norm(vecF);
			needBT = ( f >= fPrev - prm.fallTol );
		endwhile
		if ( ~doMainLoop );
			break;
		endif
		%
		%
		matInfoB(iterCount+1,:) = [ norm(vecDelta0), norm(vecDelta)/norm(vecDelta0), ...
		  norm(vecFPrev-vecF), fPrev-f, 1.0-(f/fPrev), ...
		  size(linsolfDatOut.matV,2), norm(linsolfDatOut.vecRho) ];
		iterCount++;
		matRecordX(:,iterCount+1) = vecX;
		matInfoA(iterCount+1,:) = [ iterCount, fevalCount, norm(vecX), norm(vecF) ];
		if ( f < fBest )
			fBest = f;
			vecXBest = vecX;
		endif
	endwhile
	if ( nargout>=4 )
		datOut.elapsedTime = time() - startTime;
		datOut.matRecordX = matRecordX;
		datOut.matInfoA = matInfoA;
		datOut.matInfoB = matInfoB;
	endif
return;
endfunction
