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
	fevalCount = 1;
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	matX = vecX0;
	rvecF = f0;
	matG = vecG0;
	iterCount = 0;
	%indexOfBest = 1;
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
	stepFactor = 1.0;
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
		iterCount++;
		%
		%
		if ( 1 == size(matX,2) )
			vecDelta = -stepFactor * prm.smallFallTrgt * fBest * vecGBest / sumsq(vecGBest);
		else
			%matD = matX - matX(:,indexOfBest);
			%matV = utorthdrop([ matD(:,1:indexOfBest-1), matD(:,indexOfBest+1:end) ]);
			%matX
			matD = matX - vecXBest;
			matV = utorthdrop( matD, prm.dropThresh );
			matVTD = matV' * matD;
			matVTG = matV' * matG;
			[ fFit, vecVTGFit, matVTHVFit ] = hessfit( matVTD, rvecF, matVTG );
			if ( size(matV,2) >= sizeX )
				echo__matHFit = matV * matVTHVFit * (matV')
			endif
			%
			%matB = diag(max( abs(matVTD), [], 2 ));
			vecDeltaNewton =  -stepFactor * matV * ( matVTHVFit \ vecVTGFit );
			vecGPerp = vecGBest - matV * ( matV'*vecGBest );
			vecDeltaGrad = -stepFactor * prm.gradStepCoeff * vecGPerp;
			vecDelta = vecDeltaNewton + vecGPerp;
		endif
		%
		vecXTrial = vecXBest + vecDelta;
		[ fTrial, vecGTrial ] = funchFG( vecXTrial );
		fevalCount++;
		%
		matX = [ matX, vecXTrial ];
		rvecF = [ rvecF, fTrial ];
		matG = [ matG, vecGTrial ];
		if ( fTrial < fBest )
			vecXBestPrev = vecXBest;
			fBestPrev = fBest;
			vecGBestPrev = vecGBest;
			vecXBest = vecXTrial;
			fBest = fTrial;
			vecGBest = vecGTrial;
			doStepLoop = false;
			stepFactor = 1.0;
		else
			stepFactor *= prm.btCoeff
		endif
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
return;	
endfunction


function prm = __init( funchFG, vecX0, prmIn )
	mydefs;
	sizeX = size(vecX0,1);
	prm.verbLev = VERBLEV__DETAILED;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.progressReportInterval = 0.0;
	%
	prm.fTol = eps;
	prm.gNormTol = eps;
	prm.iterLimit = 100;
	prm.fevalLimit = -1;
	prm.timeLimit = 10.0;
	prm.stopSignalCheckInterval = 1.0;
	% fall thresholds?
	%
	prm.smallFallTrgt = 0.1;
	prm.gradStepCoeff = 0.01;
	prm.dropThresh = eps^0.7;
	prm.btCoeff = 0.1;
	%
	prm = overwritefields( prm, prmIn );
	assert( prm.smallFallTrgt > 0.0 );
return;
endfunction
