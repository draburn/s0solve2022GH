function [ vecX, retCode, datOut ] = sgsolve( funchFG, init_vecX, prm=[] )
	msg( __FILE__, __LINE__, "To-do:" );
	msg( __FILE__, __LINE__, "  * Properly handle TR." );
	msg( __FILE__, __LINE__, "  * Properly set seed_vecP on limited jump. (To what?)" );
	msg( __FILE__, __LINE__, "  * Properly limit record size." );
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
	%
	seed_vecX = init_vecX;
	seed_vecP = zeros(sizeX,1);
	%
	[ sprout_vecX, sprout_vecP, sptDat ] = __evalSuperPt( funchFG, seed_vecX, seed_vecP, prm );
	fevalCount += sptDat.fevalCount;
	%
	% Set iter dat.
	seed_vecX = sprout_vecX;
	seed_vecP = sprout_vecP;
	matX = sptDat.vecXSPt;
	matG = sptDat.vecGSPt;
	rvecF = sptDat.fSPt;
	iterCount = 1;
	%
	% Set other dat.
	prev_vecX = [];
	prev_vecG = [];
	prev_f = [];
	vecX = sptDat.vecXSPt;
	vecG = sptDat.vecGSPt;
	f = sptDat.fSPt;
	%
	if ( prm.progressReportInterval >= 0.0 )
		sgsolve__reportProg;
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
		enableJumps = true;
		if (enableJumps)
		sizeR = size(matX,2);
		if ( sizeR >= 1 )
			vecXAnchor = matX(:,1);
			vecGAnchor = matG(:,1);
			fAnchor = f(:,1);
			matD = matX - vecXAnchor;
			[ matV, rvecDrop ] = utorthdrop( matD, sqrt(eps) );
			sizeK = size( matV, 2 );
			if ( sizeK >= 1 )
				%
				matVTDWC = matV' * [ matD(:,~rvecDrop), zeros(sizeX,1) ];
				matVTGWC = matV' * [ matG(:,~rvecDrop), vecGAnchor ];
				rvecFWC = [ rvecF(~rvecDrop), fAnchor ];
				fitPrm = [];
				fitPrm.epsHRegu = 0.0;
				[ fFit, vecGammaFit, matHFit ] = hessfit( matVTDWC, rvecFWC, matVTGWC, fitPrm );
				%
				doCompare = false;
				if (doCompare)
					msg( __FILE__, __LINE__, "BEGIN INFODUMP..." );
					compare_f = [ fFit, f ]
					compare_g = [ vecGammaFit, matV'*vecG ]
					matHFit
					matVTHFSV = matV'*prm.matHSecret*matV
					assert( reldiff( fFit, f ) < sqrt(eps) );
					assert( reldiff( vecGammaFit, matV'*vecG ) < sqrt(eps) );
					assert( reldiff( matHFit, matVTHFSV ) < sqrt(eps) );
					msg( __FILE__, __LINE__, "END INFODUMP." );
				endif
				%
				trSize = sqrt(sum(sum(matD.^2)));
				%trSize = 0.001*sqrt(max( sum(matD.^2,1) ));
				%trSize = sqrt(max( sum(matD.^2,1) ));
				levPrm = [];
				vecZ = levsol_eig( fFit, vecGammaFit, matHFit, [], trSize, levPrm );
				seed_vecX = vecXAnchor + matV*vecZ;
				%%%vecGammaAtZ = vecGammaFit + matHFit * vecZ; % Will be zero for full Newton step.
				vecPPerp = seed_vecP - matV*(matV'*seed_vecP);
				%%%seed_vecP = vecPPerp - matV * (vecGammaAtZ * prm.learningRate / ( 1.0 - prm.momentumFactor ));
				seed_vecP = vecPPerp;
			endif
			%
			if ( sizeK >= min([ 20, sizeX ]) )
				matX = matX( :, ~rvecDrop);
				matG = matG( :, ~rvecDrop );
				rvecF = rvecF( ~rvecDrop );
			endif
		endif
		endif
		%
		[ sprout_vecX, sprout_vecP, sptDat ] = __evalSuperPt( funchFG, seed_vecX, seed_vecP, prm );
		fevalCount += sptDat.fevalCount;
		%
		doExactnessCheck = false;
		if (doExactnessCheck)
			[ fCheck, vecGCheck ] = funchFG( sptDat.vecXSPt );
			looky_f = [ sptDat.fSPt, fCheck, sptDat.fSPt-fCheck ]
			looky_vecG = [ sptDat.vecGSPt, vecGCheck, sptDat.vecGSPt-vecGCheck ]
			error( "HALT!" );
		endif
		%
		% Move to the step.
		seed_vecX = sprout_vecX;
		seed_vecP = sprout_vecP;
		matX = [ sptDat.vecXSPt, matX ];
		matG = [ sptDat.vecGSPt, matG ];
		rvecF = [ sptDat.fSPt, rvecF ];
		prev_vecX = vecX;
		prev_vecG = vecG;
		prev_f = f;
		vecX = sptDat.vecXSPt;
		vecG = sptDat.vecGSPt;
		f = sptDat.fSPt;
		%
		% Report progress.
		if ( prm.progressReportInterval >= 0.0 && time() - progressReportedTime >= prm.progressReportInterval )
			sgsolve__reportProg;
			progressReportedTime = time();
		endif
	endwhile
	%
	if ( prm.verbLev >= VERBLEV__MAIN )
		sgsolve__reportProg;
		progressReportedTime = time();
	endif
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


function [ vecX, vecP, datOut ] = __evalSuperPt( funchFG, vecX0, vecP0, prm )
	vecX = vecX0;
	vecP = vecP0;
	datOut = [];
	datOut.fevalCount = 0;
	vecXSum = 0.0*vecX;
	vecGSum = 0.0*vecP;
	xtgSum = 0.0;
	fSum = 0.0;
	wSum = 0.0;
	for n=1:10
		[ f, vecG ] = funchFG( vecX );
		datOut.fevalCount++;
		%
		vecXSum += vecX;
		vecGSum += vecG;
		xtgSum += (vecX'*vecG);
		fSum += f;
		wSum += 1.0;
		%
		vecP = ( prm.momentumFactor * vecP ) - ( prm.learningRate * vecG );
		vecX += vecP;
	endfor
	%
	vecXAvg = (vecXSum/wSum);
	vecGAvg = (vecGSum/wSum);
	xtgAvg = (xtgSum/wSum);
	fAvg = (fSum/wSum);
	datOut.vecXSPt = vecXAvg;
	datOut.vecGSPt = vecGAvg;
	datOut.fSPt = fAvg - 0.5*( xtgAvg - vecXAvg'*vecGAvg );
return;
endfunction
