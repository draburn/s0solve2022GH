function [ vecXBest, retCode, datOut ] = sxsolve1208( funchFG, vecXInit, prm=[] )
	% Init - Universal.
	mydefs;
	startTime = time();
	%
	% Init - Default return values.
	vecXBest = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	% Init - Parse input.
	sizeX = size(vecXInit,1);
	assert( isposintscalar(sizeX) );
	assert( isrealarray(vecXInit,[sizeX,1]) );
	[ prm, fevalCount ] = __init( funchFG, vecXInit, prm );
	if ( prm.stopSignalCheckInterval >= 0.0 )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		stopSignalCheckedTime = time();
	endif
	[ fInit, vecGInit ] = funchFG( vecXInit );
	fevalCount++;
	assert( isrealscalar(fInit) );
	assert( isrealarray(vecGInit,[sizeX,1]) );
	%
	% Init - Solver.
	vecXBest = vecXInit;
	fBest = fInit;
	vecGBest = vecGInit;
	matX = [];
	rvecF = [];
	matG = [];
	iterCount = 0;
	stepCount = 0;
	currentBTFactor = 1.0;
	stepBTFactor = 0.0;
	vecXBestPrev = [];
	fBestPrev = [];
	vecGBestPrev = [];
	sizeKMostRecent = 0;
	%
	% Init - Progress reporting.
	if ( prm.progressReportInterval >= 0.0 )
		msg( __FILE__, __LINE__, sprintf( ...
		  "  %8s, %3s / %3s / %3s / %3s;  %8s, %8s;  %8s,  %8s;  %8s, %9s;  %s", ...
		  "time", ...
		  "szK", ...
		  "stp", ...
		  "rec", ...
		  "itr", ...
		  "BTFactor", ...
		  "|deltaX|", ...
		  "fBest", ...
		  "fB fall", ...
		  "|gBest|", ...
		  "|gB|fall", ...
		  "fevalCount" ) );
		msg( __FILE__, __LINE__, sprintf( ...
		  "  %8.2e, %3d / %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %9.2e;  %d", ...
		  time()-startTime, ...
		  sizeKMostRecent, ...
		  stepCount, ...
		  size(matX,2), ...
		  iterCount, ...
		  stepBTFactor, ...
		  0.0, ... % norm(vecXBestPrev - vecXBest)
		  fBest, ...
		  0.0, ... % fBestPrev - fBest
		  norm(vecGBest), ...
		  0.0, ... % norm(vecGBestPrev) - norm(vecGBest)
		  fevalCount ) );
		progressReportedTime = time();
	endif
	%
	% MAIN LOOP
	doMainLoop = true; % May be redundant.
	while (doMainLoop)
		iterCount++;
		%
		%
		% Check static stopping criteria.
		if ( prm.fTol >= 0.0 && fBest <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: f0 < prm.fTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.gTol >= 0.0 && norm(vecGBest) <= prm.gTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecGBest) <= prm.gTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.fevalLimit >= 0 && fevalCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: fevalCount >= prm.fevalLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		elseif ( prm.stepLimit >= 0 && stepCount >= prm.stepLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: stepCount >= prm.stepLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		elseif ( prm.iterLimit >= 0 && iterCount >= prm.iterLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		elseif ( prm.timeLimit >= 0.0 && time() - startTime >= prm.timeLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: time() - startTime >= prm.timeLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		elseif ( prm.stopSignalCheckInterval >= 0.0 && time() - stopSignalCheckedTime >= prm.stopSignalCheckInterval )	
			if ( stopsignalpresent() )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
				retCode = RETCODE__IMPOSED_STOP;
				doMainLoop = false;
				break;
			endif
			stopSignalCheckedTime = time();
		endif
		%
		%
		% Generate step.
		if ( 0 == size(matX,2) )
			vecDelta = -currentBTFactor * prm.gradStepCoeff * vecGBest;
			% Among other possibilities.
		else
			%[ vecDelta, datOut_getStep ] = __getStep_crude( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_simple( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_simple2( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_curateD_crude( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_curateD( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_fullSpaceFirst( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%[ vecDelta, datOut_getStep ] = __getStep_fullSpaceAlways( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			[ vecDelta, datOut_getStep ] = __getStep_bestNotFit( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			sizeKMostRecent = datOut_getStep.sizeK;
		endif
		%
		%
		% Try the step.
		vecXTrial = vecXBest + vecDelta;
		[ fTrial, vecGTrial ] = funchFG( vecXTrial );
		fevalCount++;
		%
		%
		% For step assessment, might be nice to consider median values in records?
		haveNewBest = false; % Unless...
		whollyRejectStep = false; % Unless...
		if ( fTrial > prm.horribleFCoeff * fBest )
			whollyRejectStep = true;
		elseif ( fTrial > fBest && norm(vecGTrial) > prm.horribleGCoeff * norm(vecGBest) )
			whollyRejectStep = true;
		endif
		if ( whollyRejectStep )
			currentBTFactor *= prm.btFactor;
			%msg( __FILE__, __LINE__, "Wholly rejecting step!" );
		else
			if ( fTrial < fBest )
				if (~isempty(matX))
					matD_beforeStep = matX - vecXBest; % Only for progress.
				else
					matD_beforeStep = [];
				endif
				% Put new info in *front*, so it gets orthonormalized first.
				vecXBestPrev = vecXBest;
				fBestPrev = fBest;
				vecGBestPrev = vecGBest;
				matX = [ vecXBest, matX ];
				rvecF = [ fBest, rvecF ];
				matG = [ vecGBest, matG ];
				vecXBest = vecXTrial;
				fBest = fTrial;
				vecGBest = vecGTrial;
				stepBTFactor = currentBTFactor;
				currentBTFactor = 1.0;
				haveNewBest = true;
				stepCount++;
			else
				% Put new info in *front*, so it gets orthonormalized first.
				matX = [ vecXTrial, matX ];
				rvecF = [ fTrial, rvecF ];
				matG = [ vecGTrial, matG ];
			endif
		endif
		%
		%
		% Check post-iter stop crit.
		if ( prm.deltaTol >= 0.0 && norm(vecDelta) <= prm.deltaTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecDelta) <= prm.deltaTol." );
				stepBTFactor = currentBTFactor;
				vecXBestPrev = vecXBest;
				fBestPrev = fBest;
				vecGBestPrev = vecGBest;
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		elseif ( haveNewBest )
			fFall = fBestPrev - fBest;
			if ( prm.fFallAbsTol >= 0.0 && fFall <= prm.fFallAbsTol )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: fFall <= prm.fFallAbsTol." );
				retCode = RETCODE__IMPOSED_STOP;
				doMainLoop = false;
				break;
			elseif ( prm.fFallRelTol >= 0.0 && fFall <= prm.fFallRelTol * fBestPrev )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: fFall <= prm.fFallRelTol * fBestPrev." );
				retCode = RETCODE__IMPOSED_STOP;
				doMainLoop = false;
				break;
			endif
		endif
		%
		%
		% Report progress.
		if ( haveNewBest && prm.progressReportInterval >= 0.0 )
		if ( time() - progressReportedTime >= prm.progressReportInterval )
			if (~isempty(matD_beforeStep))
				stepBeta = sqrt(max( abs(vecDelta'*matD_beforeStep)./sum( matD_beforeStep.^2, 1 ) ));
			else
				stepBeta = 1.0;
			endif
			msg( __FILE__, __LINE__, sprintf( ...
			  "  %8.2e, %3d / %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %9.2e;  %d", ...
			  time()-startTime, ...
			  sizeKMostRecent, ...
			  stepCount, ...
			  size(matX,2), ...
			  iterCount, ...
			  stepBTFactor, ...
			  norm(vecDelta), ...
			  fBest, ...
			  fBestPrev - fBest, ...
			  norm(vecGBest), ...
			  norm(vecGBestPrev) - norm(vecGBest), ...
			  fevalCount ) );
			progressReportedTime = time();
		endif
		endif
	endwhile
	%
	if ( prm.verbLev >= VERBLEV__MAIN )
		msg( __FILE__, __LINE__, sprintf( ...
		  "  %8.2e, %3d / %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %9.2e;  %d", ...
		  time()-startTime, ...
		  sizeKMostRecent, ...
		  stepCount, ...
		  size(matX,2), ...
		  iterCount, ...
		  stepBTFactor, ...
		  norm(vecDelta), ...
		  fBest, ...
		  fBestPrev - fBest, ...
		  norm(vecGBest), ...
		  norm(vecGBestPrev) - norm(vecGBest), ...
		  fevalCount ) );
		progressReportedTime = time();
	endif
	datOut.matX = matX;
	datOut.matG = matG;
	datOut.rvecF = rvecF;
	datOut.vecXBest = vecXBest;
	datOut.vecGBest = vecGBest;
	datOut.fBest = fBest;
return;	
endfunction


function [ prm, fevalCount ] = __init( funchFG, vecXInit, prmIn )
	mydefs;
	%
	% Common stuff.
	sizeX = size(vecXInit,1);
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 0.0;
	%
	% Stopping criteria - pre-step.
	prm.fTol = eps^0.5;
	prm.gTol = eps^0.5;
	prm.stepLimit = 100;
	prm.iterLimit = 100;
	prm.fevalLimit = -1;
	prm.timeLimit = 10.0;
	prm.stopSignalCheckInterval = 10.0;
	%
	% Stopping criteria - post-step.
	prm.deltaTol = eps;
	prm.fFallAbsTol = eps;
	prm.fFallRelTol = eps;
	%
	% Step generation params.
	prm.smallFallTrgt = 0.01;
	prm.gradStepCoeff = 0.01;
	%%%prm.dropThresh = eps^0.7; % Too strict!
	prm.dropThresh = eps^0.5;
	%%%prm.dropThresh = 1.0e-4;
	%prm.epsFNegativityCoeff = 0.01;
	prm.epsB = eps^0.5;
	%prm.trDCoeff = 3.0;
	prm.trDCoeff = 100.0;
	%%%prm.trDCoeff = 1e8;
	%
	% Step assessment and BT/dynamic TR params.
	prm.horribleFCoeff = 2.0;
	prm.horribleGCoeff = 2.0;
	prm.btFactor = 0.1;
	%
	% Super-point param.
	%prm.numPtPerSuperPt = 100;
	%prm.momentumCoeff = 0.9;
	%prm.learningRate = 0.01;
	%
	prm = overwritefields( prm, prmIn );
	assert( prm.smallFallTrgt > 0.0 );
	assert( prm.smallFallTrgt < 1.0 );
	assert( prm.gradStepCoeff > 0.0 );
	assert( prm.dropThresh > 0.0 );
	assert( prm.dropThresh < 1.0 );
	assert( prm.horribleGCoeff > 1.0 );
	assert( prm.horribleFCoeff > 1.0 );
	assert( prm.btFactor > 0.0 );
	assert( prm.btFactor < 1.0 );
	%
	fevalCount = 0;
return;
endfunction


function [ vecDelta, datOut ] = __getStep_crude( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	datOut.sizeK = size(matV,2);
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecZ = mycholdiv( matH, -currentBTFactor*vecGamma, false );
	vecDeltaGrad = -currentBTFactor * prm.gradStepCoeff * vecGradPerp;
	vecDelta = vecDeltaGrad + matV*vecZ;
	if (0)
		vecDScale = max( abs(matVTDWB), [], 2 );
		vecB = 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) );
		matB = diag(vecB);
		[ norm(vecZ), norm(matB*vecZ), norm(vecDelta) ]
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_simple2( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	%
	% Generate fit.
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	if ( reldiff(matV'*matV,eye(size(matV,2))) > sqrt(eps) )
		matV = utorthdrop_debug( matD, prm.dropThresh );
		echo__matD = matD
		echo__matV = matV
		echo__matVTV = matV'*matV
		error( "HALT!" );
	endif
	assert( reldiff(matV'*matV,eye(size(matV,2))) < sqrt(eps) );
	datOut.sizeK = size(matV,2);
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	%
	%
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecDeltaGradPerp = -currentBTFactor * prm.gradStepCoeff * vecGradPerp;
	%
	%
	vecDScale = max( abs(matVTDWB), [], 2 );
	matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
	bMax = currentBTFactor * prm.trDCoeff;
	
	
	%%%vecZInSpace = myhessmin( fFit, vecGamma, matH, matB, bMax );
	vecZInSpace = myhessmin( max([fFit, fBest]), vecGamma, matH, matB, bMax );
	
	
	vecDelta = vecDeltaGradPerp + matV * vecZInSpace;
	if (0)
		msg( __FILE__, __LINE__, "Infodump..." );
		[ norm(vecZInSpace), norm(matB*vecZInSpace), norm(vecDelta) ]
		echo__matB = matB
		echo__bMax = bMax
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_simple( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	%%%msg( __FILE__, __LINE__, "**********************" );
	datOut = [];
	% Generate fit.
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	if ( reldiff(matV'*matV,eye(size(matV,2))) > sqrt(eps) )
		matV = utorthdrop_debug( matD, prm.dropThresh );
		echo__matD = matD
		echo__matV = matV
		echo__matVTV = matV'*matV
		error( "HALT!" );
	endif
	assert( reldiff(matV'*matV,eye(size(matV,2))) < sqrt(eps) );
	datOut.sizeK = size(matV,2);
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	%
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecDeltaGradPerp = -currentBTFactor * prm.gradStepCoeff * vecGradPerp;
	%
	% Set scaling / (TR/)boundary matrix.
	trBeta = currentBTFactor * prm.trDCoeff;
	vecDScale = max( abs(matVTDWB), [], 2 );
	vecB = 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) );
	matB = diag(vecB);
	matBTBInv = diag(1.0./(vecB.^2));
	% We'll require ||B*z|| <= trBeta.
	%
	% Dog-leg.
	vecZSD = -( matBTBInv * vecGamma );
	zthz = vecZSD' * matH * vecZSD;
	assert( zthz > 0.0 ); % Non-positive-definite case not (yet) handled.
	%
	%%% THIS IS WRONG: vecZCauchy = (sumsq(vecZSD) / zthz) * vecZSD; <<< THAT IS WRONG.
	vecZCauchy = ( -(vecZSD'*vecGamma) / zthz) * vecZSD;
	%
	vecZNewton = mycholdiv( matH, -vecGamma, false );
	%prm_mycholdiv = [];
	%prm_mycholdiv.debugMode = true;
	%vecZNewton = mycholdiv( matH, -vecGamma, false, prm_mycholdiv );
	%
	%%%[ norm(matB*vecZCauchy), trBeta, norm(matB*vecZNewton), (matB*vecZNewton)'*(matB*vecZCauchy) ]
	%assert( norm(matB*vecZNewton) >= norm(matB*vecZCauchy)*(1.0-100.0*eps) ); % Not true if matH is perturbed to be pos-def.
	%
	if ( trBeta <= norm(matB*vecZCauchy) )
		%msg( __FILE__, __LINE__, "Taking sub-Cauchy step." );
		vecZ = vecZCauchy * trBeta / norm(matB*vecZCauchy);
	elseif ( trBeta >= norm(matB*vecZNewton) )
		%msg( __FILE__, __LINE__, "Taking full Newton step." );
		vecZ = vecZNewton;
	else
		%msg( __FILE__, __LINE__, "Taking intermediate step." );
		vecZLeg = vecZNewton - vecZCauchy;
		% a*s^2 + b*s + c = 0
		a = sumsq( matB * vecZLeg );
		b = 2.0 * (matB * vecZCauchy)' * ( matB*vecZLeg);
		c = sumsq( matB * vecZCauchy ) - trBeta^2;
		discrim = (b^2) - (4.0*a*c);
		assert( discrim >= 0.0 );
		assert( a > 0.0 );
		assert( c <= 0.0 );
		s = (-b+sqrt(discrim))/(2.0*a);
		vecZ = vecZCauchy + s*vecZLeg;
		assert( reldiff( norm(matB*vecZ) , trBeta ) < sqrt(sqrt(eps)) );
	endif
	%
	vecDelta = vecDeltaGradPerp + matV * vecZ;
	if (0)
		msg( __FILE__, __LINE__, "Infodump..." );
		[ norm(vecZ), norm(matB*vecZ), norm(vecDelta) ]
		echo__matB = matB
		echo__trBeta = trBeta
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	%
	% Generate fit.
	sizeX = size(vecXBest,1);
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	sizeK = size(matV,2);
	assert( reldiff(matV'*matV,eye(size(matV,2))) < sqrt(eps) );
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	if ( reldiff(matV'*matV,eye(size(matV,2))) > sqrt(eps) )
		echo__matD = matD
		echo__matV = matV
		echo__rd = reldiff(matV'*matV,eye(size(matV,2)))
		matV = utorthdrop_debug( matD, prm.dropThresh );
		echo__matV = matV
		echo__matVTV = matV'*matV
		echo__rd = reldiff(matV'*matV,eye(size(matV,2)))
		error( "HALT!" );
	endif
	datOut.sizeK = sizeK;
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	%
	vecDScale = max( abs(matVTDWB), [], 2 );
	matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
	bMax = currentBTFactor * prm.trDCoeff;
	%
	% Examine part of gradient that is outside of our fit space.
	vecGPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	if ( norm(vecGPerp) > prm.dropThresh * norm(vecGBest) )
		% Grab some info before we modify anything.
		h_rmsdiag = sqrt(sum(diag(matH).^2)/sizeK);
		%
		% Expand matV...
		vecVNew = vecGPerp / norm(vecGPerp);
		matV = [ matV, vecVNew ];
		%
		% Expand vecGamma...
		vecGamma = [ vecGamma; vecVNew'*vecGBest ];
		%
		% Expand non-diagonal part of matH...
		matY = matVTDWB(:,1:end-1);
		vecGommo = matG'*vecVNew - vecGamma(end);
		% matH(+) = [ matH, vecEta; vecEta, 0.0 ];
		% matY' * vecEta = vecGommo
		%  => vecEta = (matY*(matY')) \ matY * vecGommo.
		% We could do this inexactly instead.
		vecEta = mycholdiv( matY * (matY'), matY * vecGommo );
		matH = [ matH, vecEta; vecEta', 0.0 ];
		%
		% Set expanded diagonal of matH...
		% We have no information about the last element of matH.
		% But, we can at least take a few guesses...
		h_eta = sum(abs(vecEta)); % Enough to make it pos-semi-def.
		h_f = vecGPerp'*vecGPerp / (2.0*fBest); % Enough so that it doesn't go negative on its own.
		matH(end,end) = max([ h_rmsdiag, h_eta, h_f ]);
		%
		% Expand trust region...
		vecDScale = [ vecDScale; prm.gradStepCoeff / prm.trDCoeff ];
		matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
		bMax = currentBTFactor * prm.trDCoeff;
		%
		% And, go!
		vecZ = myhessmin_gperp( max([fFit, fBest]), vecGamma, matH, matB, bMax );
		vecDelta = matV * vecZ;
	else
		vecZ = myhessmin( max([fFit, fBest]), vecGamma, matH, matB, bMax );
		vecDelta = matV * vecZ;
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_curateD_crude( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	matD = matX - vecXBest;
	%
	% "Curate" matD...
	rvecDInitial = sqrt(sum(matD.^2,1));
	dTol0 = 1.0e-6;
	dTol1 = 1.0e-3;
	%
	%echo__matD = matD
	rvecKeep = logical(zeros(size(rvecF)));
	[ strongestD, strongestPt ] = max( sqrt(sum(matD.^2,1)) );
	assert( strongestD > 0.0 );
	matV = matD(:,strongestPt) / strongestD;
	rvecKeep(strongestPt) = true;
	for n = 1 : min(size(matD))
		%echo__n = n
		%echo__matV = matV
		%echo__rvecKeep = rvecKeep
		matDPerp = matD;
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		rvecDRemain = sqrt(sum(matDPerp.^2,1));
		[ strongestMerit, strongestPt ] = max( rvecDRemain./(dTol0 + (dTol1 * rvecDInitial)) );
		strongestD = rvecDRemain(strongestPt);
		%echo__rvecRemain = rvecDRemain
		%echo__norm = dTol0 + dTol1 * rvecDInitial
		%echo__merit = rvecDRemain ./  echo__norm
		%echo__strongestMerit = strongestMerit
		%echo__strongetPt = strongestPt
		if ( strongestMerit >  1.0 )
			assert( false == rvecKeep(strongestPt) );
			rvecKeep(strongestPt) = true;
			matV = [ matV, matDPerp(:,strongestPt)/strongestD ];
			if ( size(matV,2) == min(size(matD)) )
				break;
			endif
		else
			break;
		endif
	endfor
	%echo__rvecKeep = rvecKeep
	%echo__matV = matV
	%assert( max(size(matD)) < 6 );
	%
	matD = matD(:,rvecKeep);
	rvecF = rvecF(:,rvecKeep);
	matG = matG(:,rvecKeep);
	%
	datOut.sizeK = size(matV,2);
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecZ = mycholdiv( matH, -currentBTFactor*vecGamma, false );
	vecDeltaGrad = -currentBTFactor * prm.gradStepCoeff * vecGradPerp;
	vecDelta = vecDeltaGrad + matV*vecZ;
	if (0)
		vecDScale = max( abs(matVTDWB), [], 2 );
		vecB = 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) );
		matB = diag(vecB);
		[ norm(vecZ), norm(matB*vecZ), norm(vecDelta) ]
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_curateD( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	matD = matX - vecXBest;
	%
	% "Curate" matD...
	rvecDInitial = sqrt(sum(matD.^2,1));
	dTol0 = 1.0e-6;
	dTol1 = 1.0e-3;
	%
	%echo__matD = matD
	rvecKeep = logical(zeros(size(rvecF)));
	[ strongestD, strongestPt ] = max( sqrt(sum(matD.^2,1)) );
	assert( strongestD > 0.0 );
	matV = matD(:,strongestPt) / strongestD;
	rvecKeep(strongestPt) = true;
	for n = 1 : min(size(matD))
		%echo__n = n
		%echo__matV = matV
		%echo__rvecKeep = rvecKeep
		matDPerp = matD;
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		rvecDRemain = sqrt(sum(matDPerp.^2,1));
		[ strongestMerit, strongestPt ] = max( rvecDRemain./(dTol0 + (dTol1 * rvecDInitial)) );
		strongestD = rvecDRemain(strongestPt);
		%echo__rvecRemain = rvecDRemain
		%echo__norm = dTol0 + dTol1 * rvecDInitial
		%echo__merit = rvecDRemain ./  echo__norm
		%echo__strongestMerit = strongestMerit
		%echo__strongetPt = strongestPt
		if ( strongestMerit >  1.0 )
			assert( false == rvecKeep(strongestPt) );
			rvecKeep(strongestPt) = true;
			matV = [ matV, matDPerp(:,strongestPt)/strongestD ];
			if ( size(matV,2) == min(size(matD)) )
				break;
			endif
		else
			break;
		endif
	endfor
	%echo__rvecKeep = rvecKeep
	%echo__matV = matV
	%assert( max(size(matD)) < 6 );
	%
	matD = matD(:,rvecKeep);
	rvecF = rvecF(:,rvecKeep);
	matG = matG(:,rvecKeep);
	%
	sizeK = size(matV,2);
	datOut.sizeK = sizeK;
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	%
	%if ( sizeK == size(vecXBest,1) )
	%	echo__matHSecret = prm.matHSecret
	%	echo__matVHVT = matV*matH*matV'
	%	%error("HALT!");
	%endif
	%
	vecDScale = max( abs(matVTDWB), [], 2 );
	matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
	bMax = currentBTFactor * prm.trDCoeff;
	%
	% Examine part of gradient that is outside of our fit space.
	vecGPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	if ( norm(vecGPerp) > prm.dropThresh * norm(vecGBest) )
		assert( size(matV,2) < size(matV,1) );
		% Grab some info before we modify anything.
		h_rmsdiag = sqrt(sum(diag(matH).^2)/sizeK);
		%
		% Expand matV...
		vecVNew = vecGPerp / norm(vecGPerp);
		matV = [ matV, vecVNew ];
		%
		% Expand vecGamma...
		vecGamma = [ vecGamma; vecVNew'*vecGBest ];
		%
		% Expand non-diagonal part of matH...
		matY = matVTDWB(:,1:end-1);
		vecGommo = matG'*vecVNew - vecGamma(end);
		% matH(+) = [ matH, vecEta; vecEta, 0.0 ];
		% matY' * vecEta = vecGommo
		%  => vecEta = (matY*(matY')) \ matY * vecGommo.
		% We could do this inexactly instead.
		vecEta = mycholdiv( matY * (matY'), matY * vecGommo );
		matH = [ matH, vecEta; vecEta', 0.0 ];
		%
		% Set expanded diagonal of matH...
		% We have no information about the last element of matH.
		% But, we can at least take a few guesses...
		h_eta = sum(abs(vecEta)); % Enough to make it pos-semi-def.
		h_f = vecGPerp'*vecGPerp / (2.0*fBest); % Enough so that it doesn't go negative on its own.
		matH(end,end) = max([ h_rmsdiag, h_eta, h_f ]);
		%
		% Expand trust region...
		vecDScale = [ vecDScale; prm.gradStepCoeff / prm.trDCoeff ];
		matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
		bMax = currentBTFactor * prm.trDCoeff;
		%
		% And, go!
		vecZ = myhessmin_gperp( max([fFit, fBest]), vecGamma, matH, matB, bMax );
		vecDelta = matV * vecZ;
	else
		vecZ = myhessmin( max([fFit, fBest]), vecGamma, matH, matB, bMax );
		vecDelta = matV * vecZ;
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_fullSpaceFirst( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	sizeX = size(matX,1);
	numRecords = size(matX,2);
	if ( numRecords < sizeX )
		assert( currentBTFactor > 0.99 );
		vecDelta = zeros(sizeX,1);
		vecDelta(numRecords+1) = -0.1 * currentBTFactor * prm.gradStepCoeff * vecGBest(numRecords+1);
		datOut.sizeK = 0;
		return;
	endif
	if ( numRecords == sizeX )
		vecXAnchor = vecXBest;
		vecGAnchor = vecGBest;
		fAnchor = fBest;
	else
		vecXAnchor = matX(:,numRecords-sizeX);
		vecGAnchor = matG(:,numRecords-sizeX);
		fAnchor = rvecF(numRecords-sizeX);
	endif
	matX = matX(:,1+numRecords-sizeX:end);
	matG = matG(:,1+numRecords-sizeX:end);
	rvecF = rvecF(1+numRecords-sizeX:end);
	%
	matD = matX - vecXAnchor;
	matV = orth( matD );
	matVTDWB = matV'*[ matD, zeros(sizeX,1) ];
	matVTGWB = matV'*[ matG, vecGAnchor ];
	rvecFWB = [ rvecF, fAnchor ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	%vecDelta = -currentBTFactor * matV * ( matH \ vecGamma );
	rd = reldiff( matV*matH*(matV'), prm.matHSecret );
	vecDelta = matV * mycholdiv( matH, -currentBTFactor*(matV'*vecGBest) );
	datOut.sizeK = size(matV,2);
	
	if (0)
		%echo__vecXBest = vecXBest
		%echo__matX = matX
		%echo__matD = matD
		%echo__matV = matV
		%echo__matVTDWB = matVTDWB
		%cond(matVTDWB(:,1:end-1))
		%echo__matVHVT = matV*matH*(matV')
		echo__matHRes = matV*matH*(matV') - prm.matHSecret
		%echo__vecXNext = vecXBest + vecDelta
		error( "HALT!" );
	endif
return;
endfunction

function [ vecDelta, datOut ] = __getStep_fullSpaceAlways( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	sizeX = size(matX,1);
	%
	matX = [];
	matG = [];
	rvecF = [];
	for n = 1 : sizeX
		v = zeros(sizeX,1);
		v(n) = 1.0;
		vecX = vecXBest + v;
		[ f, vecG ] = prm.funchFGSecret( vecX );
		rvecF(n) = f;
		matX(:,n) = vecX;
		matG(:,n) = vecG;
	end
	matD = matX - vecXBest;
	matV = orth( matD );
	matVTDWB = matV'*[ matD, zeros(sizeX,1) ];
	matVTGWB = matV'*[ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGammaFit, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecDelta = matV * mycholdiv( matH, -currentBTFactor*(matV'*vecGBest) );
	datOut.sizeK = size(matV,2);
	if (0)
		%echo__vecXBest = vecXBest
		%echo__matD = matD
		%echo__matV = matV
		%echo__matVTDWB = matVTDWB
		%cond(matVTDWB(:,1:end-1))
		%echo__matVHVT = matV*matH*(matV')
		echo__matHRes = matV*matH*(matV') - prm.matHSecret
		%echo__vecXNext = vecXBest + vecDelta
		error( "HALT!" );
	endif
return;
endfunction
function [ vecDelta, datOut ] = __getStep_bestNotFit( currentBTFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	matD = matX - vecXBest;
	%
	% "Curate" matD...
	rvecDInitial = sqrt(sum(matD.^2,1));
	dTol0 = 1.0e-6;
	dTol1 = 1.0e-3;
	%
	%echo__matD = matD
	rvecKeep = logical(zeros(size(rvecF)));
	[ strongestD, strongestPt ] = max( sqrt(sum(matD.^2,1)) );
	assert( strongestD > 0.0 );
	matV = matD(:,strongestPt) / strongestD;
	rvecKeep(strongestPt) = true;
	for n = 1 : min(size(matD))
		%echo__n = n
		%echo__matV = matV
		%echo__rvecKeep = rvecKeep
		matDPerp = matD;
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		matDPerp = matDPerp - (matV * ( matV'*matDPerp ));
		rvecDRemain = sqrt(sum(matDPerp.^2,1));
		[ strongestMerit, strongestPt ] = max( rvecDRemain./(dTol0 + (dTol1 * rvecDInitial)) );
		strongestD = rvecDRemain(strongestPt);
		%echo__rvecRemain = rvecDRemain
		%echo__norm = dTol0 + dTol1 * rvecDInitial
		%echo__merit = rvecDRemain ./ (dTol0 + dTol1 * rvecDInitial)
		%echo__rvecKeep = rvecKeep
		%echo__strongestMerit = strongestMerit
		%echo__strongetPt = strongestPt
		if ( strongestMerit >  1.0 )
			assert( false == rvecKeep(strongestPt) );
			rvecKeep(strongestPt) = true;
			matV = [ matV, matDPerp(:,strongestPt)/strongestD ];
			if ( size(matV,2) == min(size(matD)) )
				break;
			endif
		else
			break;
		endif
	endfor
	%echo__rvecKeep = rvecKeep
	%echo__matV = matV
	%assert( max(size(matD)) < 6 );
	%
	matD = matD(:,rvecKeep);
	rvecF = rvecF(:,rvecKeep);
	matG = matG(:,rvecKeep);
	%
	sizeK = size(matV,2);
	datOut.sizeK = sizeK;
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGammaFit, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecGammaBest = matV'*vecGBest;
	%
	%if ( sizeK == size(vecXBest,1) )
	%	echo__matHSecret = prm.matHSecret
	%	echo__matVHVT = matV*matH*matV'
	%	%error("HALT!");
	%endif
	%
	vecDScale = max( abs(matVTDWB), [], 2 );
	matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
	bMax = currentBTFactor * prm.trDCoeff;
	%
	% Examine part of gradient that is outside of our fit space.
	vecGPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	if ( norm(vecGPerp) > prm.dropThresh * norm(vecGBest) )
		assert( size(matV,2) < size(matV,1) );
		% Grab some info before we modify anything.
		h_rmsdiag = sqrt(sum(diag(matH).^2)/sizeK);
		%
		% Expand matV...
		vecVNew = vecGPerp / norm(vecGPerp);
		matV = [ matV, vecVNew ];
		%
		% Expand vecGammaBest...
		vecGammaBest = [ vecGammaBest; vecVNew'*vecGBest ];
		%
		% Expand non-diagonal part of matH...
		matY = matVTDWB(:,1:end-1);
		vecGommo = matG'*vecVNew - vecGammaBest(end);
		% matH(+) = [ matH, vecEta; vecEta, 0.0 ];
		% matY' * vecEta = vecGommo
		%  => vecEta = (matY*(matY')) \ matY * vecGommo.
		% We could do this inexactly instead.
		vecEta = mycholdiv( matY * (matY'), matY * vecGommo );
		matH = [ matH, vecEta; vecEta', 0.0 ];
		%
		% Set expanded diagonal of matH...
		% We have no information about the last element of matH.
		% But, we can at least take a few guesses...
		h_eta = sum(abs(vecEta)); % Enough to make it pos-semi-def.
		h_f = vecGPerp'*vecGPerp / (2.0*fBest); % Enough so that it doesn't go negative on its own.
		matH(end,end) = max([ h_rmsdiag, h_eta, h_f ]);
		%
		% Expand trust region...
		vecDScale = [ vecDScale; prm.gradStepCoeff / prm.trDCoeff ];
		matB = diag( 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) ) );
		bMax = currentBTFactor * prm.trDCoeff;
		%
		% And, go!
		vecZ = myhessmin_gperp( max([fFit, fBest]), vecGammaBest, matH, matB, bMax );
		vecDelta = matV * vecZ;
	else
		vecZ = myhessmin( max([fFit, fBest]), vecGammaBest, matH, matB, bMax );
		vecDelta = matV * vecZ;
	endif
return;
endfunction