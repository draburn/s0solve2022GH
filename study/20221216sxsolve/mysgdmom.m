function [ vecX, retCode, datOut ] = mysgdmom( funchFG, init_vecX, prm=[] )
	% Init - Universal.
	mydefs;
	startTime = time();
	%
	% Init - Default return values.
	vecX = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	% Init - Parse input.
	sizeX = size(init_vecX,1);
	assert( isposintscalar(sizeX) );
	assert( isrealarray(init_vecX,[sizeX,1]) );
	[ prm, fevalCount ] = __init( funchFG, init_vecX, prm );
	if ( prm.stopSignalCheckInterval >= 0.0 )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		stopSignalCheckedTime = time();
	endif
	[ init_f, init_vecG ] = funchFG( init_vecX );
	fevalCount++;
	assert( isrealscalar(init_f) );
	assert( isrealarray(init_vecG,[sizeX,1]) );
	%
	% Init - Solver.
	vecX = init_vecX;
	vecG = init_vecG;
	f = init_f;
	vecDelta = zeros(sizeX,1);
	iterCount = 0;
	%
	% Init - Extras.
	prev_vecX = [];
	prev_vecG = [];
	prev_f = [];
	%
	if ( prm.progressReportInterval >= 0.0 )
		mysgdmom__reportProg;
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
		if ( prm.fTol >= 0.0 && f <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: f < prm.fTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.gTol >= 0.0 && norm(vecG) <= prm.gTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecG) <= prm.gTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.xTol >= 0.0 && ~isempty(prev_vecX) && norm(vecX-prev_vecX) <= prm.xTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecX-prev_vecX) <= prm.xTol." );
			retCode = RETCODE__IMPOSED_STOP;
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
		vecDelta = ( prm.momentumFactor * vecDelta ) - ( prm.learningRate * vecG );
		%
		%
		% Try the step.
		trial_vecX = vecX + vecDelta;
		[ trial_f, trial_vecG ] = funchFG( trial_vecX );
		fevalCount++;
		assert( isrealscalar(trial_f) );
		assert( isrealarray(trial_vecG,[sizeX,1]) );
		%
		% Move to the step.
		prev_vecX = vecX;
		prev_vecG = vecG;
		prev_f = f;
		vecX = trial_vecX;
		vecG = trial_vecG;
		f = trial_f;
		%
		% Report progress.
		if ( prm.progressReportInterval >= 0.0 && time() - progressReportedTime >= prm.progressReportInterval )
			mysgdmom__reportProg;
			progressReportedTime = time();
		endif
	endwhile
	%
	if ( prm.verbLev >= VERBLEV__MAIN )
		mysgdmom__reportProg;
		progressReportedTime = time();
	endif
	datOut.vecG = vecG;
	datOut.f = f;
	datOut.prm = prm;
return;	
endfunction


function [ prm, fevalCount ] = __init( funchFG, vecXInit, prmIn )
	mydefs;
	%
	% Common stuff.
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 1.0;
	%
	% Stopping criteria.
	prm.xTol = eps^0.7;
	prm.fTol = eps^0.7;
	prm.gTol = eps^0.4;
	prm.iterLimit = 1E6;
	prm.fevalLimit = 1E6;
	prm.timeLimit = 30.0;
	prm.stopSignalCheckInterval = 3.0;
	%
	% General step generation param.
	prm.learningRate = 0.1;
	prm.momentumFactor = 0.9;
	%
	prm = overwritefields( prm, prmIn );
	%
	fevalCount = 0;
return;
endfunction
