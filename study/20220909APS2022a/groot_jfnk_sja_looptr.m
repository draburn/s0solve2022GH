% jfnk_sja_looptr = _baseline + sja + in-(linsolf)-loop trust region

function [ vecXBest, grootFlag, fevalCount, datOut ] = groot_jfnk_sja_looptr( funchF, vecX0, prm=[] )
	
	msg( __FILE__, __LINE__, "WARNING: WORK-IN-PROGRESS!!!" );
	
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
	if ( prm.verbLev >= VERBLEV__DETAILS )
		msg( __FILE__, __LINE__, "Welcome!" );
	endif
	%
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, "Updating of approximate Jacobian on step is not implemented!" );
	%
	vecF0 = funchF( vecX0 );
	fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	f0 = norm(vecF0);
	%
	%
	useStepSearch = mygetfield( prm, "useStepSearch", true );
	assert( isbool(useStepSearch) );
	assert( isscalar(useStepSearch) );
	if ( useStepSearch )
		stepType = mygetfield( prm, "stepType", "powell" );
		btCoeff = mygetfield( prm, "btCoeff", 0.2 );
		assert( isrealscalar(btCoeff) );
		assert( 0.0 < btCoeff );
		useTR = mygetfield( prm, "useTR", true );
		assert( isbool(useTR) );
		assert( isscalar(useTR) );
		if ( useTR )
			ftCoeff = mygetfield( prm, "ftCoeff", 1.5 );
			assert( isrealscalar(ftCoeff) );
			assert( 0.0 < ftCoeff );
		endif
	endif
	useAP = mygetfield( prm, "useAP", true );
	assert( isbool(useAP) );
	assert( isscalar(useAP) );
	if ( useAP )
		apUpdateType = mygetfield( prm, "apUpdateType", "secant" );
		useWoodbury = mygetfield( prm, "useWoodbury", true );
		assert( isbool(useWoodbury) );
		assert( isscalar(useWoodbury) );
	endif
	useSJA = mygetfield( prm, "useSJA", true );
	assert( isbool(useSJA) );
	assert( isscalar(useSJA) );
	if ( useSJA )
		sja_matJAInv = [];
	endif
	%
	%
	datOut.prm = prm;
	vecX = vecX0;
	vecF = vecF0;
	f = f0;
	stepInverseTRLimit = 0.0; % Limit step (||delta||) to Trust Region size.
	iterCount = 0;
	matRecordX(:,iterCount+1) = vecX;
	matInfoA(iterCount+1,:) = [ iterCount, fevalCount, norm(vecX), norm(vecF) ];
	matInfoB = [];
	fBest = f0;
	vecXBest = vecX0;
	if (useAP)
	if (~useWoodbury)
		assert( sizeF == sizeX );
		matJA = eye(sizeX);
	else
		assert( sizeF == sizeX );
		matJAInv = eye(sizeX);
	endif
	endif
	doMainLoop = true;
	while (doMainLoop)
		if ( prm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( "  %3d,  %6d:  %10.3e.", iterCount, fevalCount, f ) );
		endif
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
		%
		% NEW FOR LOOPTR
		assert( useStepSearch );
		assert( useTR );
		linsolfPrm.btCoeff = mygetfield( prm, "btCoeff", btCoeff );
		linsolfPrm.stepType = stepType;
		linsolfPrm.stepInverseTRLimit = stepInverseTRLimit;
		%
	if ( useSJA && ~isempty(sja_matJAInv) )
		linsolfPrm.matP = sja_matJAInv;
		linsolfPrm.useSJA = false;
		[ vecDelta0, linsolfDatOut ] = linsolf_sja_looptr( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
	else
		if (useAP)
		if (~useWoodbury)
			linsolfPrm.matP = pinv(matJA); % Ugly.
			[ vecDelta0, linsolfDatOut ] = linsolf_sja_looptr( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		else
			linsolfPrm.matP = matJAInv;
			[ vecDelta0, linsolfDatOut ] = linsolf_sja_looptr( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		endif
		else
			[ vecDelta0, linsolfDatOut ] = linsolf_sja_looptr( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		endif
	endif
		fevalCount += linsolfDatOut.fevalCount;
		if ( fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		endif
		if ( ~isrealarray(vecDelta0,[sizeX,1]) )
			msgif( prm.verbLev >= VERBLEV__FLAG, __FILE__, __LINE__, "INTERNAL ERROR: linsolf_sja_looptr() returned an invalid value." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		endif
		sizeK = size(linsolfDatOut.matV,2);
		%
		%
		vecXPrev = vecX;
		vecFPrev = vecF;
		fPrev = f;
		%
		%
		vecDelta = vecDelta0;
		doSearchLoop = true;
		if ( ~useStepSearch || ~useTR )
			% If not using TR, then reset every time.
			stepInverseTRLimit = 0.0;  % Possibly unnecessary code, to be safe.
		endif
		%
		%
		% NEW FOR LOOPTR
		vecX = vecXPrev + vecDelta;
		vecF = funchF( vecX );
		fevalCount++;
		f = norm(vecF);
		%
		%
		if ( useStepSearch && useTR )
			stepInverseTRLimit = 1.0/( ftCoeff * norm(vecDelta) );
		else
			stepInverseTRLimit = 0; % Possibly unnecessary code, to be safe.
		endif
		assert( norm(vecDelta) > 0.0 );
		if ( useAP )
		if (~useWoodbury)
			switch ( tolower(apUpdateType) )
			case { "none" }
				matJA += ( linsolfDatOut.matW - (matJA*linsolfDatOut.matV) ) * (linsolfDatOut.matV');
			case { "full space secant" }
				matJA += ( linsolfDatOut.matW - (matJA*linsolfDatOut.matV) ) * (linsolfDatOut.matV');
				matJA += ( vecF - vecFPrev - (matJA*vecDelta) ) * ((vecDelta')/sumsq(vecDelta));
			case { "secant", "broyden" }
				matW = linsolfDatOut.matW;
				matV = linsolfDatOut.matV;
				vecY = matV'*vecDelta;
				matW += ( vecF - vecFPrev - matW*vecY ) * ((vecY')/sumsq(vecDelta));
				matJA += ( matW - (matJA*matV) ) * (matV');
			case { "osqu" }
				% "On Step Quadratic Update"
				matW = linsolfDatOut.matW;
				matV = linsolfDatOut.matV;
				vecY = matV'*vecDelta;
				matW += 2.0*( vecF - vecFPrev - matW*vecY ) * ((vecY')/sumsq(vecDelta));
				matJA += ( matW - (matJA*matV) ) * (matV');
			otherwise
				error([ "Invalid value of apUpdateType (\"" apUpdateType "\"." ]);
			endswitch
		else
			matV = linsolfDatOut.matV;
			matW = linsolfDatOut.matW;
			vecY = matV'*vecDelta;
			switch ( tolower(apUpdateType) )
			case { "none" }
				matWUpdated = matW;
				%
				%matD = (matJAInv * matWUpdated) - matV;
				%matA = eye(sizeX) + matD * (matV');
				%matAInv = eye(sizeX) - matD * pinv( eye(sizeK) + matV'*matD ) * (matV');
				%assert( reldiff( matA * matAInv, eye(sizeX) ) < sizeX*sqrt(eps) );
				%assert( reldiff( matAInv * matA, eye(sizeX) ) < sizeX*sqrt(eps) );
				%matJAInv = matAInv * matJAInv;
			case { "secant", "broyden" }
				matWUpdated = matW + ( vecF - vecFPrev - matW*vecY ) * ((vecY')/sumsq(vecDelta));
			case { "osqu" } % "On Step Quadratic Update"
				matWUpdated = matW + 2.0*( vecF - vecFPrev - matW*vecY ) * ((vecY')/sumsq(vecDelta));
			otherwise
				error([ "Invalid value of apUpdateType (\"" apUpdateType "\"." ]);
			endswitch
			matD = (matJAInv * matWUpdated) - matV;
			matIK = eye(sizeK);
			matJAInv -= ( matD * ( (matIK + matV'*matD) \ (matV') ) ) * matJAInv;
		endif
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
		if ( useSJA && isempty(sja_matJAInv) )
			sja_matJAInv = linsolfDatOut.sja_matJAInv;
		endif
	endwhile
	if ( nargout>=4 )
		datOut.elapsedTime = time() - startTime;
		datOut.matRecordX = matRecordX;
		datOut.matInfoA = matInfoA;
		datOut.matInfoB = matInfoB;
	endif
	if ( prm.verbLev >= VERBLEV__DETAILS )
		msg( __FILE__, __LINE__, "Goodbye!" );
	endif
return;
endfunction
