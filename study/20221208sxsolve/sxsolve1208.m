function [ vecXBest, retCode, datOut ] = sxsolve1208( funchFG, vecX0, prm=[] )
	mydefs;
	startTime = time();
	vecXBest = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	sizeX = size(vecX0,1);
	assert( isposintscalar(sizeX) );
	assert( isrealarray(vecX0,[sizeX,1]) );
	prm = __init( funchFG, vecX0, prm );
	if ( prm.stopSignalCheckInterval >= 0.0 )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		stopSignalCheckedTime = time();
	endif
	%
	[ f0, vecG0 ] = funchFG( vecX0 );
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	matX = vecX0;
	rvecF = f0;
	matG = vecG0;
	fevalCount = 1; % Includes geval.
	iterCount = 0;
	%
	vecXBestPrev = vecX0;
	fBestPrev = f0;
	vecGBestPrev = vecG0;
	vecXBest = vecX0;
	fBest = f0;
	vecGBest = vecG0;
	%
	progressReportedTime = time();
	doMainLoop = true;
	while (doMainLoop)
		if ( prm.fTol >= 0.0 && fBest <= prm.fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: f0 < prm.fTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.gNormTol >= 0.0 && norm(vecGBest) <= prm.gNormTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: norm(vecGBest) <= prm.gNormTol." );
			retCode = RETCODE__SUCCESS;
			doMainLoop = false;
			break;
		elseif ( prm.fevalLimit >= 0 && iterCount >= prm.fevalLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.fevalLimi." );
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
		if ( prm.progressReportInterval >= 0.0 )
		if ( 0 == iterCount || time() - progressReportedTime >= prm.progressReportInterval )
			msg( __FILE__, __LINE__, sprintf( ...
			  " %8.2e, %3d;  %2d;  %8.2e;  %8.2e,  %8.2e;  %8.2e, %8.2e; %d", ...
			  time()-startTime, ...
			  iterCount, ...
			  size(matX,2), ...
			  norm(vecXBestPrev - vecXBest), ...
			  fBest, ...
			  fBestPrev - fBest, ...
			  norm(vecGBest), ...
			  norm(vecGBestPrev) - norm(vecGBest), ...
			  fevalCount ) );
			progressReportedTime = time();
		endif
		endif
		iterCount++
		%
		%
		% DO WORK HERE.
		%
		%
	endwhile
	%
	if ( prm.verbLev >= VERBLEV__MAIN )
		msg( __FILE__, __LINE__, sprintf( ...
		  " %8.2e, %3d;  %2d;  %8.2e;  %8.2e,  %8.2e;  %8.2e, %8.2e; %d", ...
		  time()-startTime, ...
		  iterCount, ...
		  size(matX,2), ...
		  norm(vecXBestPrev - vecXBest), ...
		  fBest, ...
		  fBestPrev - fBest, ...
		  norm(vecGBest), ...
		  norm(vecGBestPrev) - norm(vecGBest), ...
		  fevalCount ) );
		progressReportedTime = time();
	endif
	%
return;	
endfunction


function prm = __init( funchFG, vecX0, prmIn )
	mydefs;
	sizeX = size(vecX0,1);
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 3.0;
	%
	prm.fTol = eps;
	prm.gNormTol = eps;
	prm.iterLimit = 1;
	prm.fevalLimit = -1;
	prm.timeLimit = 10.0;
	prm.stopSignalCheckInterval = 1.0;
	% fall thresholds?
	%
	prm = overwritefields( prm, prmIn );
return;
endfunction