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
	currentTRFactor = 1.0;
	vecXBestPrev = [];
	fBestPrev = [];
	vecGBestPrev = [];
	%
	% Init - Progress reporting.
	if ( prm.progressReportInterval >= 0.0 )
		msg( __FILE__, __LINE__, sprintf( ...
		  "  %8s, %3s / %3s / %3s;  %8s, %8s;  %8s,  %8s;  %8s, %8s;  %s", ...
		  "time", ...
		  "stp", ...
		  "szX", ...
		  "itr", ...
		  "TRFactor", ...
		  "|deltaX|", ...
		  "fBest", ...
		  "fB fall", ...
		  "|gBest|", ...
		  "|gB|fall", ...
		  "fevalCount" ) );
		msg( __FILE__, __LINE__, sprintf( ...
		  "  %8.2e, %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %8.2e;  %d", ...
		  time()-startTime, ...
		  stepCount, ...
		  size(matX,2), ...
		  iterCount, ...
		  currentTRFactor, ...
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
			vecDelta = -currentTRFactor * prm.gradStepCoeff * vecGBest;
			% Among other possibilities.
		else
			vecDelta = __getStep_crude( currentTRFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
			%vecDelta = __getStep( currentTRFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
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
			currentTRFactor *= prm.btFactor;
		else
			if ( fTrial < fBest )
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
				currentTRFactor = 1.0;
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
			msg( __FILE__, __LINE__, sprintf( ...
			  "  %8.2e, %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %8.2e;  %d", ...
			  time()-startTime, ...
			  stepCount, ...
			  size(matX,2), ...
			  iterCount, ...
			  currentTRFactor, ...
			  norm(vecXBestPrev - vecXBest), ...
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
		  "  %8.2e, %3d / %3d / %3d;  %8.2e, %8.2e;  %8.2e,  %8.2e;  %8.2e, %8.2e;  %d", ...
		  time()-startTime, ...
		  stepCount, ...
		  size(matX,2), ...
		  iterCount, ...
		  currentTRFactor, ...
		  norm(vecXBestPrev - vecXBest), ...
		  fBest, ...
		  fBestPrev - fBest, ...
		  norm(vecGBest), ...
		  norm(vecGBestPrev) - norm(vecGBest), ...
		  fevalCount ) );
		progressReportedTime = time();
	endif
return;	
endfunction


function [ prm, fevalCount ] = __init( funchFG, vecXInit, prmIn )
	mydefs;
	%
	% Common stuff.
	sizeX = size(vecXInit,1);
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 1.0;
	%
	% Stopping criteria - pre-step.
	prm.fTol = eps;
	prm.gTol = eps;
	prm.stepLimit = 999;
	prm.iterLimit = 999;
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
	prm.smallFallTrgt = 0.1;
	prm.gradStepCoeff = 0.1;
	prm.dropThresh = eps^0.7;
	%prm.epsFNegativityCoeff = 0.01;
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


function [ vecDelta, datOut ] = __getStep_crude( currentTRFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecDeltaNewton = matV * mycholdiv( matH, -currentTRFactor*vecGamma );
	vecDeltaGrad = -currentTRFactor * prm.gradStepCoeff * vecGradPerp;
	vecDelta = vecDeltaNewton + vecDeltaGrad;
	datOut = [];
return;
endfunction
function [ vecDelta, datOut ] = __getStep( currentTRFactor, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	% Generate fit.
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	
	error( "Perhaps should do step without adding gradPerp first?" );
	
	
	% HAXOR
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecDeltaNewton = matV * mycholdiv( matH, -currentTRFactor*vecGamma );
	vecDeltaGrad = -currentTRFactor * prm.gradStepCoeff * vecGradPerp;
	vecDelta = vecDeltaNewton + vecDeltaGrad;
	%msg( __FILE__, __LINE__, "..." );
	%matG
	%vecG_model = matV * ( vecGamma + matH * (matV'*matD) )
	if ( false && size(matV,2) >= size(matV,1)-1 )
		msg( __FILE__, __LINE__, "..." );
		matHFS = matV * matH * (matV');
		resH = sqrt( ...
		  sum(sum((matHFS-prm.matHSecret).^2)) ...
		 / ( sum(sum(matHFS.^2)) + sum(sum(prm.matHSecret.^2)) ) )
		matH_fullspace = matV * matH * (matV')
		prm.matHSecret
	endif
	
	%
	%
	% (Possibly) add one more vector to V: the gradient direction from the launch point.
	% This version of the code (2022-12-09) assumes the launch point is always the best point,
	%  as well as the "anchor" for the subspace basis.
	% But, there's an additional assumption here: we *could* do some fit over our records,
	%  but, instead, we'll simply the gradient calculated at this point.
	vecT = vecGBest;
	vecTPerp = vecT - ( matV * ( matV' * vecT ) );
	if ( norm(vecTPerp) > prm.dropThresh * norm(vecT) )
		matY = matVTDWB(:,1:end-1);
		vecVNew = vecTPerp / norm(vecTPerp);
		%
		matV = [ matV, vecVNew ];
		vecGamma = [ vecGamma; vecVNew'*vecGBest ];
		vecGommo = matG'*vecVNew - vecGamma(end);
		% matH(+) = [ matH, vecEta; vecEta, 0.0 ];
		% matY' * vecEta = vecGommo
		%  => vecEta = (matY*(matY')) \ matY * vecGommo.
		% TODO: We could do this inexactly instead.
		vecEta = mycholdiv( matY * (matY'), matY * vecGommo );
		matH = [ matH, vecEta; vecEta', 0.0 ];
		
		%msg( __FILE__, __LINE__, "..." );
		%vecG_model = matV * ( vecGamma + matH * (matV'*matD) )
		if ( false && size(matV,2) >= size(matV,1)-1 )
			msg( __FILE__, __LINE__, "..." );
			matHFS = matV * matH * (matV');
			resH = sqrt( ...
			  sum(sum((matHFS-prm.matHSecret).^2)) ...
			 / ( sum(sum(matHFS.^2)) + sum(sum(prm.matHSecret.^2)) ) )
			matH_fullspace = matV * matH * (matV')
			prm.matHSecret
		endif
	else
		% Too small. Can't add it.
	endif
	return
	%
	% This code is crude.
	vecDeltaNewton = matV * mycholdiv( matH, -currentTRFactor*vecGamma );
	%vecDeltaGrad = -currentTRFactor * prm.gradStepCoeff * vecGradPerp;
	%vecDelta = vecDeltaNewton + vecDeltaGrad;
	
	
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecDeltaNewton = matV * mycholdiv( matH, -currentTRFactor*vecGamma );
	%vecDeltaGrad = -currentTRFactor * prm.gradStepCoeff * vecGradPerp;
	%vecDelta = vecDeltaNewton + vecDeltaGrad;
	datOut = [];
	%
return;
endfunction
