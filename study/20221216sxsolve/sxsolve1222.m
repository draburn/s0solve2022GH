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
	trSize = [];
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
		[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		sizeK = datOut.sizeK;
		%
		%
		% Try the step.
		vecXTrial = vecXBest + vecDelta;
		[ fTrial, vecGTrial ] = funchFG( vecXTrial );
		fevalCount++;
		%
		%
		% For step assessment, might be nice to consider median values in records?
		trSize_beforeUpdate = trSize;
		if ( fTrial > prm.horribleCoeff * fBest )
			trialCount_horrible++;
			haveNewBest = false;
			trSize = min([ trSize, norm(vecDelta) ]) * prm.btFactor_horrible;
		elseif ( fTrial >= fBest )
			trialCount_bad++;
			% Put new info in *front*, so it gets orthonormalized first.
			matX = [ vecXTrial, matX ];
			rvecF = [ fTrial, rvecF ];
			matG = [ vecGTrial, matG ];
			haveNewBest = false;
			trSize = min([ trSize, norm(vecDelta) ]) * prm.btFactor_bad;
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
			haveNewBest = true;
			trSize = max([ trSize, norm(vecDelta)*prm.ftFactor_good ]);
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
		if ( haveNewBest && prm.progressReportInterval >= 0.0 && time() - progressReportedTime >= prm.progressReportInterval )
			sxsolve1222__reportProg;
			progressReportedTime = time();
		endif
		clear haveNewBest;
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
	prm.progressReportInterval = 0.01;
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
	prm.defaultFFallTarget = 0.01;
	prm.horribleCoeff = 2.0;
	prm.btFactor_horrible = 0.1;
	prm.btFactor_bad = 0.5;
	prm.ftFactor_good = 1.3;
	prm.deltaTol = eps^0.8;
	%
	prm = overwritefields( prm, prmIn );
	%
	fevalCount = 0;
return;
endfunction


function [ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	%defaultGradCoeff = 0.1;
	%vecDelta = -defaultGradCoeff*vecGBest;
	%
	if ( isempty(trSize) )
		defaultFallTarget = 0.01;
		vecDelta = (-defaultFallTarget*fBest/sumsq(vecGBest)) * vecGBest;
	else
		vecDelta = (-trSize/norm(vecGBest)) * vecGBest;
	endif
	datOut.sizeK = 0;
return;
endfunction
