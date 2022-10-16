% jfnk_sja = _baseline + sja

function [ vecXBest, grootFlag, fevalCount, datOut ] = groot_jfnk_sja( funchF, vecX0, prm=[] )
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
	if ( useSJA && ~isempty(sja_matJAInv) )
		linsolfPrm.matP = sja_matJAInv;
		linsolfPrm.useSJA = false;
		[ vecDelta0, linsolfDatOut ] = linsolf_sja( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
	else
		if (useAP)
		if (~useWoodbury)
			linsolfPrm.matP = pinv(matJA); % Ugly.
			[ vecDelta0, linsolfDatOut ] = linsolf_sja( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		else
			linsolfPrm.matP = matJAInv;
			[ vecDelta0, linsolfDatOut ] = linsolf_sja( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		endif
		else
			[ vecDelta0, linsolfDatOut ] = linsolf_sja( funchMatJProd, -vecF, zeros(sizeX,1), linsolfPrm );
		endif
	endif
		fevalCount += linsolfDatOut.fevalCount;
		if ( fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
			grootFlag = GROOT_FLAG__STOP;
			doMainLoop = false; % Redundant.
			break;
		endif
		sizeK = size(linsolfDatOut.matV,2);
		% Subspace Newton step...
		vecYN = linsolfDatOut.vecY;
		% Subspace Cauchy step...
		% Take y = -p*matS*matW'*vecF, for scaling matrix matS.
		% Find value of p such that d/dp ( 0.5*||vecF_model||^2 ) = 0.
		matS = eye(sizeK);
		vecG = linsolfDatOut.matW'*vecF;
		assert( norm(vecG) > 0.0 );
		vecGS =  matS*vecG;
		assert( norm(vecGS) > 0.0 );
		p = vecGS'*vecG / sumsq( linsolfDatOut.matW * vecGS );
		vecYC = -p*vecGS;
		assert( norm(vecYC) > 0.0 );
		vecYD = vecYN - vecYC;
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
		while (doSearchLoop)
			msgif( prm.verbLev >= VERBLEV__DETAILED, __FILE__, __LINE__, sprintf( ...
			  "  BT  %3d,  %6d:  %10.3e,  %10.3e,  %10.3e.", ...
			  iterCount, fevalCount, f, norm(vecDelta), stepInverseTRLimit ) );
			if ( norm(vecDelta) * stepInverseTRLimit <= 1.0 + 100.0*eps )
				% Step is in TR.
				vecX = vecXPrev + vecDelta;
				vecF = funchF( vecX );
				fevalCount++;
				f = norm(vecF);
				if ( f <= fPrev - prm.fallTol )
					doSearchLoop = false; % Possibly unnecessary code, to be safe.
					break;
				endif
				if ( ~useStepSearch )
					msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Reached fallTol and doStepSearch is false." );
					grootFlag = GROOT_FLAG__FAIL;
					doMainLoop = false;
					doSearchLoop = false; % Possibly unnecessary code, to be safe.
					break;
				endif
				targetStepSize = btCoeff * norm(vecDelta);
			else
				targetStepSize = ( 1.0 - 100.0*eps ) / stepInverseTRLimit;
			endif
			%
			assert( useStepSearch );
			if ( norm(vecDelta) <= prm.stepTol )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Reached stepTol." );
				grootFlag = GROOT_FLAG__FAIL;
				doMainLoop = false;
				doSearchLoop = false; % Possibly unnecessary code, to be safe.
				break;
			elseif ( fevalCount >= prm.fevalLimit )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached fevalLimit." );
				grootFlag = GROOT_FLAG__STOP;
				doMainLoop = false;
				doSearchLoop = false; % Possibly unnecessary code, to be safe.
				break;
			elseif ( stopsignalpresent() )
				msgif( prm.verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
				grootFlag = GROOT_FLAG__STOP;
				doMainLoop = false;
				doSearchLoop = false; % Possibly unnecessary code, to be safe.
				break;
			endif
			%
			switch ( tolower(stepType) )
			case { "newton", "newt" }
				vecDelta = targetStepSize * vecDelta0 / norm(vecDelta0);
			case { "powell", "dogleg", "dog leg", "dog-leg" }
				if ( targetStepSize < norm(vecYC) )
					vecDelta = targetStepSize * linsolfDatOut.matV * ( vecYC / norm(vecYC) );
				else
					% Find where the second leg intersects the boundary.
					a = sumsq(vecYD);
					b = 2.0*(vecYC'*vecYD);
					c = sumsq(vecYC) - targetStepSize^2;
					discrim = (b^2) - (4.0*a*c);
					assert( discrim >= 0.0 );
					assert( 0.0 < a );
					p = (-b+sqrt(discrim))/(2.0*a);
					vecDelta = linsolfDatOut.matV * ( vecYC + p*vecYD );
					assert( abs(norm(vecDelta)-targetStepSize) < 0.01*targetStepSize );
				endif
			otherwise
				error([ "Invalid value of stepType (\"" stepType "\")." ]);
			endswitch
		endwhile
		if ( ~doMainLoop );
			break;
		endif
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
