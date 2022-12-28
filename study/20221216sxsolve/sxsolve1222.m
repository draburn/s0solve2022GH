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
	trSize = prm.initialTRSize; % Could be empty.
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
		%[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		%[ vecDelta, datOut ] = __getStep_newtCheat( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		%[ vecDelta, datOut ] = __getStep_crude( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		%[ vecDelta, datOut ] = __getStep_crude_north( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		if (1)
		if ( 0 == mod(iterCount,2) )
			[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		else
			[ vecDelta, datOut ] = __getStep_crude_north( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		endif
		endif
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
		if (0) %%%%%%%%%%%%%%%%%%%%
		if ( fTrial > prm.horribleCoeff * fBest )
			trialCount_horrible++;
			haveNewBest = false;
			if ( isempty(trSize) )
				trSize = norm(vecDelta) * prm.btFactor_horrible;
			else
				trSize *= prm.btFactor_horrible;
			endif
		elseif ( fTrial >= fBest )
			trialCount_bad++;
			% Put new info in *front*, so it gets orthonormalized first.
			matX = [ vecXTrial, matX ];
			rvecF = [ fTrial, rvecF ];
			matG = [ vecGTrial, matG ];
			haveNewBest = false;
			if ( isempty(trSize) )
				trSize = norm(vecDelta) * prm.btFactor_bad;
			else
				trSize *= prm.btFactor_bad;
			endif
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
			if (~isempty(trSize))
				if ( 0 == trialCount_bad && 0 == trialCount_horrible )
					trSize = max([ trSize, norm(vecDelta)*prm.ftFactor_goodStart ]);
				else
					trSize = max([ trSize, norm(vecDelta)*prm.ftFactor_good ]);
				endif
			endif
		endif
		else %%%%%%%%%%%%%%%%%%%%
		if ( fTrial > prm.horribleCoeff * fBest )
			trialCount_horrible++;
			if ( isempty(trSize) )
				trSize = norm(vecDelta) * prm.btFactor_horrible;
			else
				trSize *= prm.btFactor_horrible;
			endif
		elseif ( fTrial >= fBest )
			trialCount_bad++;
			if ( isempty(trSize) )
				trSize = norm(vecDelta) * prm.btFactor_bad;
			else
				trSize *= prm.btFactor_bad;
			endif
		else
			trialCount_good++;
			if (~isempty(trSize))
				if ( 0 == trialCount_bad && 0 == trialCount_horrible )
					trSize = max([ trSize, norm(vecDelta)*prm.ftFactor_goodStart ]);
				else
					trSize = max([ trSize, norm(vecDelta)*prm.ftFactor_good ]);
				endif
			endif
		endif
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
		endif %%%%%%%%%%%%%%%%%%%%
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
	prm.progressReportInterval = 0.0;
	%
	% Stopping criteria - pre-step.
	prm.fTol = eps^0.7;
	prm.gTol = eps^0.4;
	prm.iterLimit = 1E4;
	prm.fevalLimit = 80;
	prm.timeLimit = 1000.0;
	prm.stopSignalCheckInterval = 3.0;
	%
	% General step generation param.
	prm.maxFallTarget = 0.1;
	%prm.initialTRSize = 1.0e-3;
	prm.initialTRSize = [];
	prm.horribleCoeff = 2.0;
	prm.btFactor_horrible = 0.1;
	prm.btFactor_bad = 0.5; % Could even be 1.0, assuming matX is used.
	prm.ftFactor_good = 1.5;
	prm.ftFactor_goodStart = 10.0;
	prm.deltaTol = eps^0.8;
	%
	prm = overwritefields( prm, prmIn );
	%
	fevalCount = 0;
return;
endfunction


function [ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	vecDelta = (-prm.maxFallTarget*fBest/sumsq(vecGBest)) * vecGBest;
	if ( ~isempty(trSize) && norm(vecDelta) > trSize )
		vecDelta *= trSize / norm(vecDelta);
	endif
	datOut.sizeK = 0;
return;
endfunction
function [ vecDelta, datOut ] = __getStep_newtCheat( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	datOut = [];
	sizeX = size(vecXBest,1);
	matV = eye(sizeX,sizeX);
	matH = zeros(sizeX,sizeX);
	epsFD = 1.0e-4;
	for n=1:sizeX
		xp = vecXBest + epsFD*matV(:,n);
		xm = vecXBest - epsFD*matV(:,n);
		[ fp, gp ] = prm.funchFGSecret( xp );
		[ fm, gm ] = prm.funchFGSecret( xm );
		matH(:,n) = (gp-gm)/(2.0*epsFD);
	endfor
	matH = ( matH' + matH ) / 2.0;
	vecDelta = newtish_eig( matH, -vecGBest );
	if ( ~isempty(trSize) && norm(vecDelta) > trSize )
		vecDelta *= trSize / norm(vecDelta);
	endif
	datOut.sizeK = sizeX;
return;
endfunction
function [ vecDelta, datOut ] = __getStep_crude( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	vecDelta = [];
	datOut = [];
	sizeX = size(vecXBest,1);
	%
	if ( 0 == size(matX,2) )
		[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		return;
	endif
	%
	matD = matX - vecXBest;
	[ matV, rvecDrop ] = utorthdrop( matD, 0.1 );
	sizeK = size(matV,2);
	%
	vecGammaBest = matV'*vecGBest;
	vecGPerp = vecGBest - matV*vecGammaBest;
	if ( norm(vecGPerp) > 0.5 * norm(vecGBest) )
		[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		return;
	endif
	%
	matVTDWB = [ matV'*matD(:,~rvecDrop), zeros(sizeK,1) ];
	matVTGWB = [ matV'*matG(:,~rvecDrop), vecGammaBest ];
	rvecFWB = [ rvecF(~rvecDrop), fBest ];
	fitPrm = [];
	fitPrm.epsHRegu = 0.0;
	[ fFit, vecGammaFit, matHFit ] = hessfit( matVTDWB, rvecFWB, matVTGWB, fitPrm );
	%
	validateFit = true;
	if ( validateFit )
		%lambdaVec = eig(matHFit);
		%eigRangeHFit = [ min(lambdaVec), max(lambdaVec) ]
		matHTrue = prm.matHSecret;
		rdF = reldiff(fFit,fBest);
		rdG = reldiff(vecGammaFit,vecGammaBest);
		rdH = reldiff(matHFit,matV'*matHTrue*matV);
		if ( rdF >= 1.0e-4 || rdG >= 1.0e-4 || rdH >= 1.0e-4 )
			%vecXBest
			%matX
			%matD
			%vecGBest
			%matG
			matV
			matVTDWB
			log10(sum(matVTDWB.^2,1))
			matVTGWB
			log10(sum(matVTGWB.^2,1))
			rvecFWB
			fVals = [ fBest, fFit, fFit - fBest ]
			vecGammaVals = [ vecGammaBest, vecGammaFit, vecGammaFit - vecGammaBest ]
			matHFit
			matVTHVTrue = matV'*matHTrue*matV
			matHRes = matHFit - matVTHVTrue
			rd = [ rdF, rdG, rdH ]
		endif
		assert( rdF < 1.0e-4 );
		assert( rdG < 1.0e-4 );
		assert( rdH < 1.0e-4 );
	endif
	vecDelta = matV * levsol_eig( fBest, vecGammaBest, matHFit, [], trSize );
	%
	datOut.sizeK = sizeK;
return;
endfunction
function [ vecDelta, datOut ] = __getStep_crude_north( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm )
	vecDelta = [];
	datOut = [];
	sizeX = size(vecXBest,1);
	%
	if ( 0 == size(matX,2) )
		[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		return;
	endif
	matD_raw = matX - vecXBest;
	matV_perpCheck = utorthdrop( matD_raw, 0.1 );
	vecGPerp = vecGBest - matV_perpCheck*(matV_perpCheck'*vecGBest);
	if ( norm(vecGPerp) > 0.1 * norm(vecGBest) )
		[ vecDelta, datOut ] = __getStep_grad( trSize, vecXBest, fBest, vecGBest, matX, rvecF, matG, prm );
		return;
	endif
	[ foo, rvecDrop ] = utorthdrop( matD_raw, sqrt(eps) );
	%%%[ foo, rvecDrop ] = utorthdrop( matD_raw, 0.001 );
	%%%rvecDrop = (1:size(matX,2))>sizeX;
	matX = matX(:,~rvecDrop);
	matG = matG(:,~rvecDrop);
	rvecF = rvecF(~rvecDrop);
	%
	matD = matX - vecXBest;
	sizeD = size(matD,2);
	matE = [ eye(sizeD), zeros(sizeD,1) ];
	vecGammaBest = matD'*vecGBest;
	matDTGWB = matD'*[ matG, vecGBest ];
	rvecFWB = [ rvecF, fBest ];
	fitPrm = [];
	fitPrm.epsHRegu = 0.0;
	[ fFit, vecGammaFit, matHFit ] = hessfit( matE, rvecFWB, matDTGWB, fitPrm );
	%
	validateFit = false;
	if ( validateFit )
		%lambdaVec = eig(matHFit);
		%eigRangeHFit = [ min(lambdaVec), max(lambdaVec) ]
		matHTrue = prm.matHSecret;
		rdF = reldiff(fFit,fBest);
		rdG = reldiff(vecGammaFit,vecGammaBest);
		rdH = reldiff(matHFit,matD'*matHTrue*matD);
		if ( rdF >= 1.0e-4 || rdG >= 1.0e-4 || rdH >= 1.0e-4 )
			fVals = [ fBest, fFit, fFit - fBest ]
			vecGammaVals = [ vecGammaBest, vecGammaFit, vecGammaFit - vecGammaBest ]
			matHFit
			matVTHVTrue = matD'*matHTrue*matD
			matHRes = matHFit - matVTHVTrue
			rd = [ rdF, rdG, rdH ]
		endif
		assert( rdF < 1.0e-4 );
		assert( rdG < 1.0e-4 );
		assert( rdH < 1.0e-4 );
	endif
	%
	levPrm = [];
	levPrm.epsEig = eps;
	matB = diag(sum(matD.^2,1));
	vecZ = levsol_eig( fFit, vecGammaFit, matHFit, matB, trSize, levPrm );
	vecDelta =  matD * vecZ;
	%
	compareFitToNoiseless = false;
	if ( compareFitToNoiseless )
		vecGTrue = prm.funch_vecGSecret( vecXBest );
		vecGTrueHat = vecGTrue / norm(vecGTrue);
		vecGBestHat = vecGBest / norm(vecGBest);
		%g_trueXBest = vecGTrueHat.'*vecGBestHat
		%
		vecGammaTrue = matD'*vecGTrue;
		vecGammaTrueHat = vecGammaTrue / norm(vecGammaTrue);
		vecGammaBestHat = vecGammaBest / norm(vecGammaBest);
		vecGammaFitHat = vecGammaFit / norm(vecGammaFit);
		%
		%gammma_trueXFit = vecGammaTrueHat'*vecGammaFitHat
		%gammma_trueXBest = vecGammaTrueHat'*vecGammaBestHat
		%gamma_fitXBest = vecGammaFitHat'*vecGammaBestHat
		%
		[ fTrue, vecGTrue ] = prm.funchFGNoiselessSecret(vecXBest);
		matHFSTrue = prm.matHSecret;
		%
		compare_f = [ fBest, fFit, fTrue, fFit-fTrue ]
		%
		vecGammaTrue = matD'*vecGTrue;
		compare_gamma = [ vecGammaBest, vecGammaFit, vecGammaTrue, vecGammaFit-vecGammaTrue ]
		rdG = reldiff(vecGammaFit,vecGammaTrue)
		%
		matHFit
		matHTrue = matD'*matHFSTrue*matD
		matHRes = matHFit - matHTrue
		rdH = reldiff(matHFit,matHTrue)
		%
		stepSize = norm(vecDelta)
		vecXNext = vecXBest + vecDelta;
		fNextFit = fFit + vecZ'*vecGammaFit + (vecZ'*matHFit*vecZ)/2.0;
		[ fNextTrue, vecGNextTrue ] = prm.funchFGNoiselessSecret(vecXNext);
		compare_fNext = [ fNextFit, fNextTrue, fNextFit-fNextTrue ]
		%gNextFit = vecGBest + ???
	endif
	%
	datOut.sizeK = size(matD,2);
endfunction
