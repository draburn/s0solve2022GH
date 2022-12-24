function [ vecXBest, retCode, datOut ] = sxsolve1222( funchFG, vecXInit, prm=[] )
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
	currentBTFactor = 1.0;
	stepBTFactor = 0.0;
	vecXBestPrev = [];
	fBestPrev = [];
	vecGBestPrev = [];
	sizeKMostRecent = 0;
	%
	% Init - Progress reporting.
	if ( prm.progressReportInterval >= 0.0 )
		sxsolve1222__reportProg;
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
		gradStepCoeff = mygetfield( prm, "gradStepCoeff", 0.1 );
		vecDelta = -currentBTFactor *gradStepCoeff * vecGBest;
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
			sxsolve1222__reportProg;
			progressReportedTime = time();
		endif
		endif
	endwhile
	%
	if ( prm.verbLev >= VERBLEV__MAIN )
		sxsolve1222__reportProg;
		progressReportedTime = time();
	endif
	datOut.matX = matX;
	datOut.matG = matG;
	datOut.rvecF = rvecF;
	datOut.vecXBest = vecXBest;
	datOut.vecGBest = vecGBest;
	datOut.fBest = fBest;
	datOut.prm = prm;
return;	
endfunction


function [ prm, fevalCount ] = __init( funchFG, vecXInit, prmIn )
	mydefs;
	%
	% Common stuff.
	sizeX = size(vecXInit,1);
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 0.1;
	%
	% Stopping criteria - pre-step.
	prm.fTol = eps^0.8;
	prm.gTol = eps^0.8;
	prm.iterLimit = 1E4;
	prm.fevalLimit = -1;
	prm.timeLimit = 30.0;
	prm.stopSignalCheckInterval = 10.0;
	%
	% Stopping criteria - post-step.
	prm.deltaTol = eps^0.8;
	%
	% Step assessment and BT/dynamic TR params.
	prm.horribleFCoeff = 2.0;
	prm.horribleGCoeff = 2.0;
	prm.btFactor = 0.1;
	%
	prm = overwritefields( prm, prmIn );
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
	datOut.sizeK = 0;
	if (0)
		vecDScale = max( abs(matVTDWB), [], 2 );
		vecB = 1.0 ./ ( vecDScale + prm.epsB * max(vecDScale) );
		matB = diag(vecB);
		[ norm(vecZ), norm(matB*vecZ), norm(vecDelta) ]
	endif
return;
endfunction
