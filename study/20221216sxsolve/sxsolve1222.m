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
	stepSizeCoeff = 1.0;
	%
	% Init - Extras.
	trialCount_horrible = 0;
	trialCount_bad = 0;
	trialCount_good = 0;
	prev_vecXBest = [];
	prev_fBest = [];
	prev_vecGBest = [];
	%
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
		vecDelta = -stepSizeCoeff * gradStepCoeff * vecGBest;
		clear gradStepCoeff;
		sizeK = 0;
		%
		%
		% Try the step.
		vecXTrial = vecXBest + vecDelta;
		[ fTrial, vecGTrial ] = funchFG( vecXTrial );
		fevalCount++;
		%
		%
		% For step assessment, might be nice to consider median values in records?
		stepSizeCoeff_beforeUpdate = stepSizeCoeff;
		if ( fTrial > prm.horribleCoeff * fBest )
			trialCount_horrible++;
			stepSizeCoeff *= prm.btFactor;
		elseif ( fTrial >= fBest )
			trialCount_bad++;
			% Put new info in *front*, so it gets orthonormalized first.
			matX = [ vecXTrial, matX ];
			rvecF = [ fTrial, rvecF ];
			matG = [ vecGTrial, matG ];
			% Should we reset stepSizeCoeff?
		else
			trialCount_good++;
			% Put new info in *front*, so it gets orthonormalized first.
			prev_vecXBest = vecXBest;
			prev_fBest = fBest;
			prev_vecGBest = vecGBest;
			matX = [ vecXBest, matX ];
			rvecF = [ fBest, rvecF ];
			matG = [ vecGBest, matG ];
			vecXBest = vecXTrial;
			fBest = fTrial;
			vecGBest = vecGTrial;
			stepSizeCoeff = 1.0;
		endif
		%
		%
		% Check post-iter stop crit.
		if ( prm.deltaTol >= 0.0 && norm(vecDelta) <= prm.deltaTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecDelta) <= prm.deltaTol." );
			prev_vecXBest = vecXBest;
			prev_fBest = fBest;
			prev_vecGBest = vecGBest;
			retCode = RETCODE__IMPOSED_STOP;
			doMainLoop = false;
			break;
		endif
		%
		%
		% Report progress.
		if ( prm.progressReportInterval >= 0.0 && time() - progressReportedTime >= prm.progressReportInterval )
			sxsolve1222__reportProg;
			progressReportedTime = time();
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
	% General step generation param.
	%%%prm.stepType = "grad";
	prm.horribleCoeff = 2.0;
	prm.btFactor = 0.1;
	prm.deltaTol = eps^0.8;
	%
	prm = overwritefields( prm, prmIn );
	%
	fevalCount = 0;
return;
endfunction


function [ vecDelta, datOut ] = __getStep_crude( stepSizeCoeff, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	matD = matX - vecXBest;
	matV = utorthdrop( matD, prm.dropThresh );
	datOut.sizeK = size(matV,2);
	matVTDWB = matV' * [ matD, zeros(size(vecXBest)) ];
	matVTGWB = matV' * [ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	[ fFit, vecGamma, matH ] = hessfit( matVTDWB, rvecFWB, matVTGWB );
	vecGradPerp = vecGBest - ( matV * ( matV' * vecGBest ) );
	vecZ = mycholdiv( matH, -stepSizeCoeff*vecGamma, false );
	vecDeltaGrad = -stepSizeCoeff * prm.gradStepCoeff * vecGradPerp;
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
