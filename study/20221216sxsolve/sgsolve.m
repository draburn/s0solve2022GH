function [ vecX, retCode, datOut ] = sgsolve( funchFG, init_vecX, prm=[] )
	if (0)
	% DRaburn 2023-01-04:
	%  This list was last update... 2023-01-02, perhaps?
	msg( __FILE__, __LINE__, "Need:" );
	msg( __FILE__, __LINE__, "  ( Testing? Timing? Refactor? )" );
	msg( __FILE__, __LINE__, "  * Meander/SPt: Reduce 'learningRate' if clearly overshooting." );
	msg( __FILE__, __LINE__, "  * Core: Fix 'best' vs 'latest'!" );
	msg( __FILE__, __LINE__, "Want:" );
	msg( __FILE__, __LINE__, "  * Core: Reduce repeated manipulation of record data (matD, matV, matVT*)." );
	msg( __FILE__, __LINE__, "  * Output: Report ledger and subspace size." );
	msg( __FILE__, __LINE__, "Maybe Consider:" );
	msg( __FILE__, __LINE__, "  * Core: Only accept jump if implies sufficient reduction to ||g||?." );
	msg( __FILE__, __LINE__, "  * Refine hessfit: From ledger dat, two- or three-pass, allow nonorthog basis." );
	msg( __FILE__, __LINE__, "  * Refine jump: Incl variation of (any/all) out-of-space components of gradient." );
	msg( __FILE__, __LINE__, "  * Refine jump: Analyze 'alpha' (handling out-of-space quants)." );
	msg( __FILE__, __LINE__, "  * Core: Curation record data (matX, matG; optim num fevals vs own work)." );
	msg( __FILE__, __LINE__, "      Consider: discard if contribution to reducing ||g|| is small." );
	msg( __FILE__, __LINE__, "  * Test trFactor (which is prop matD)." );
	msg( __FILE__, __LINE__, "  * BT/TR: Improve good/bad check on fSPt and dynamic step size limit in addition to 'trFactor' (which is prop matD)." );
	msg( __FILE__, __LINE__, "Potential Optimizations:" );
	msg( __FILE__, __LINE__, "  * Meander/SPt: Improved 'momentum' (via estimating gradient at lead point?)." );
	msg( __FILE__, __LINE__, "  * Meander/SPt: Grab data from just start & end (elim per-feval work?)." );
	endif
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
	rvecW = sptDat.wSPt;
	iterCount = 1;
	stepSizeCoeff = 1.0;
	simple_trSize = [];
	%
	% Set other dat.
	prev_vecX = [];
	prev_vecG = [];
	prev_f = [];
	init_g = sptDat.vecGSPt;
	init_f = sptDat.fSPt;
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
		if (0)
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
				%
				if (0)
					% Trial code for just pointint at Newton point. Didn't work.
					vecPInPlane = vecXAnchor + matV*vecZ - seed_vecX; % Point to in-plane newon pt.
					vecPInPlane *= norm(matV'*seed_vecP)/norm(vecPInPlane); % Scale per in-plane amount.
					seed_vecP = vecPPerp + vecPInPlane;
				endif
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
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_basicCts( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_simpleFitCts( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_simpleFitReduced( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_simple_redux( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, simple_trSize, prm );
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_simple_ineffective( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, simple_trSize, prm );
		%[ trial_vecX, trial_vecP, jumpDat ] = __jump_january( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
		[ trial_vecX, trial_vecP, jumpDat ] = __jump_januarySeedy( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
		assert( isrealarray(trial_vecX,[sizeX,1]) );
		assert( isrealarray(trial_vecP,[sizeX,1]) );
		%
		if (1)
			if (~isempty(jumpDat.rvecUseForFit))
				% 2022-12-29-1930: Crude curation criteria.
				% Note that we could alternatively *combine* points.
				matX = matX(:,jumpDat.rvecUseForFit);
				matG = matG(:,jumpDat.rvecUseForFit);
				rvecF = rvecF(jumpDat.rvecUseForFit);
				rvecW = rvecW(jumpDat.rvecUseForFit);
			endif
			ledgerLimit = mygetfield( prm, "ledgerLimit", [] );
			if ( ~isempty(ledgerLimit) && size(matX,2) > ledgerLimit )
				matX = matX(:,1:ledgerLimit);
				matG = matG(:,1:ledgerLimit);
				rvecF = rvecF(1:ledgerLimit);
				rvecW = rvecW(1:ledgerLimit);
			endif
		endif
		%
		[ sprout_vecX, sprout_vecP, sptDat ] = __evalSuperPt( funchFG, trial_vecX, trial_vecP, prm );
		fevalCount += sptDat.fevalCount;
		assert( sptDat.fSPt < init_f / eps );
		%
		doExactnessCheck = false;
		if (doExactnessCheck)
			[ fCheck, vecGCheck ] = funchFG( sptDat.vecXSPt );
			looky_f = [ sptDat.fSPt, fCheck, sptDat.fSPt-fCheck ]
			looky_vecG = [ sptDat.vecGSPt, vecGCheck, sptDat.vecGSPt-vecGCheck ]
			error( "HALT!" );
		endif
		%
		if ( sptDat.fSPt >= f )
		%%%if ( norm(sptDat.vecGSPt) >= norm(vecG) )
		%%%if (0)
			% That didn't work so well.
			simple_trSize = 0.1*norm( trial_vecX - seed_vecX );
			% Ignore the result and just do another superPt.
			%msg( __FILE__, __LINE__, sprintf( "Ouch: %g >= %g.", sptDat.fSPt, f ) );
			[ sprout_vecX, sprout_vecP, sptDat ] = __evalSuperPt( funchFG, seed_vecX, seed_vecP, prm );
			fevalCount += sptDat.fevalCount;
			assert( sptDat.fSPt < init_f / eps );
			%msg( __FILE__, __LINE__, sprintf( "Instead: %g vs %g.", sptDat.fSPt, f ) );
		else
			if ( ~isempty(simple_trSize) )
				simple_trSize = 2.0*norm( trial_vecX - seed_vecX );
			endif
			%msg( __FILE__, __LINE__, sprintf( "Yays: %g < %g.", sptDat.fSPt, f ) );
		endif
		%
		% Move to the step.
		seed_vecX = sprout_vecX;
		seed_vecP = sprout_vecP;
		matX = [ sptDat.vecXSPt, matX ];
		matG = [ sptDat.vecGSPt, matG ];
		rvecF = [ sptDat.fSPt, rvecF ];
		rvecW = [ sptDat.wSPt, rvecW];
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
	sizeX = size(vecX0,1);
	numFevalPerSuperPt = mygetfield( prm, "numFevalPerSuperPt", 10 );
	for n=1:numFevalPerSuperPt
		[ f, vecG ] = funchFG( vecX );
		assert( isrealscalar(f) );
		assert( isrealarray(vecG,[sizeX,1]) );
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
	datOut.wSPt = wSum;
return;
endfunction


function [ vecXNew, vecPNew, jumpDat ] = __jump_basicCts( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	% DRaburn 2023-01-03:
	%  Original jump routine for this solver.
	%  The "Cts" means that, in the limit of stepSizeCoeff -> 0,
	%   the output "new" values are the input "seed" values.
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	%
	% Generate subspace basis matrix.
	basisDropThresh = mygetfield( prm, "basisDropThresh", 0.01 );
	matD = matX - vecXAnchor;
	matV = utorthdrop( matD, basisDropThresh );
	jumpDat.sizeK = size(matV,2);
	if ( size(matV,2) == 0 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Drop anything that's too far outside of subspace.
	rvecDSq = sum( matD.^2, 1 );
	matVTD = matV'*matD;
	rvecVTDSq = sum( matVTD.^2, 1 );
	fitDropThresh = mygetfield( prm, "fitDropThresh", 0.001 );
	% basisDropThresh is threshold for dropping, but fitDropThresh is threshold for not dropping???
	rvecUseForFit = ( rvecVTDSq >= (1.0-fitDropThresh) * rvecDSq );
	rvecUseForFit(indexAnchor) = true;
	jumpDat.rvecUseForFit = rvecUseForFit;
	matDForFit = matD(:,rvecUseForFit);
	matGForFit = matG(:,rvecUseForFit);
	rvecFForFit = rvecF(rvecUseForFit);
	rvecWForFit = rvecW(rvecUseForFit);
	%
	% Generate fit.
	matVTDForFit = matV'*matDForFit;
	matVTGForFit = matV'*matGForFit;
	fitPrm = [];
	fitPrm.epsHRegu = 0.0;
	%msg( __FILE__, __LINE__, sprintf( "rcond(mtm(matVTDForFit)) = %g", rcond(mtm(matVTDForFit)) ) );
	%msg( __FILE__, __LINE__, sprintf( "rcond(mtm(matVTDForFit')) = %g", rcond(mtm(matVTDForFit')) ) );
	[ fFit, vecGammaFit, matHFit ] = hessfit( matVTDForFit, rvecFForFit, matVTGForFit, fitPrm );
	if ( fFit < 0.0 )
		warning( "fFit is negative." );
	endif
	
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matV'*vecGTrue;
		matHFit
		matVTHFSV = matV'*prm.matHSecret*matV
		matRes = matHFit - matVTHFSV
		rdH = reldiff( matHFit, matVTHFSV )
		vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
		rdG = reldiff( vecGammaFit, vecGammaTrue )
		fTrue
		fFit
		rdF = reldiff( fFit, fTrue )
	endif
	
	
	%
	if ( mygetfield( prm, "validateFit", false ) );
		matVTHFSV = matV'*prm.matHSecret*matV;
		rdF = reldiff( fFit, fAnchor );
		rdG = reldiff( vecGammaFit, matV'*vecGAnchor );
		rdH = reldiff( matHFit, matVTHFSV );
		rdTol = 100.0*sqrt(eps);
		if ( rdF > rdTol || rdG > rdTol || rdH > rdTol )
			msg( __FILE__, __LINE__, "BEGIN INFODUMP..." );
			compare_f = [ fFit, fAnchor, fFit-fAnchor ]
			compare_g = [ vecGammaFit, matV'*vecGAnchor, vecGammaFit-matV'*vecGAnchor ]
			matHFit
			matVTHFSV
			matHRes = matHFit - matVTHFSV
			rdVals = [ rdF, rdG, rdH ]
			msg( __FILE__, __LINE__, "END INFODUMP." );
		endif
		assert( rdF <= rdTol );
		assert( rdG <= rdTol );
		assert( rdH <= rdTol );
	endif
	%
	% Decompose (pre-jump) "seed" x.
	vecDSeed = vecXSeed - vecXAnchor;
	vecYSeed = matV'*vecDSeed;
	vecXPerp = vecDSeed - matV*vecYSeed;
	assert( reldiff( vecXSeed, vecXAnchor + (matV*vecYSeed) + vecXPerp ) <= sqrt(eps) );
	assert( norm(matV'*vecXPerp) <= sqrt(eps)*(norm(vecXAnchor)+norm(vecXSeed)) );
	%
	% These "seed" values are just estimates at the (subspace-projected) seed.
	fSeed = fFit + vecYSeed'*vecGammaFit + (vecYSeed'*matHFit*vecYSeed)/2.0;
	vecGammaSeed = vecGammaFit + matHFit*vecYSeed;
	normGammaSeed = norm(vecGammaSeed);
	%
	% Decompose (pre-jump) "seed" momentum.
	% Note that we actually do know how the full gradient varies throughout the subspace,
	%  not just how "gamma" (the part of the gradient that's in the subspace) varies;
	%  appropriate handling would essentially allow us to incorporate "vecPPerp" in to
	%  "coeffPGamma" and "vecGammaPerp".
	% Meh.
	vecT = matV'*vecPSeed;
	vecPPerp = vecPSeed - (matV*vecT);
	if ( 0.0 < normGammaSeed )
		coeffPGamma = (vecGammaSeed'*vecT) / (normGammaSeed^2);
	else
		coeffPGamma = 0.0;
		msg( __FILE__, __LINE__, "WARNING: 0.0 >= normGammaSeed; this should never happen." );
	endif
	vecGammaPerp = vecT - coeffPGamma*vecGammaSeed;
	assert( reldiff( matV*( coeffPGamma*vecGammaSeed + vecGammaPerp ) + vecPPerp, vecPSeed ) <= sqrt(eps) );
	assert( norm(matV'*vecPPerp) <= sqrt(eps)*norm(vecPSeed) );
	assert( abs(vecGammaPerp'*vecGammaSeed) <= sqrt(eps)*norm(vecPSeed) );
	%
	% Find point on Levenberg curve subject to trust region.
	%
	%% 2022-12-29-1642: Current code is a placeholder which at least obeys "stepSizeCoeff" in some form.
	%matB = [];
	%bMax = [];
	%levPrm = [];
	%% Note: We want to start the curve from the (subspace-projected) "seed", not the anchor.
	%vecZFull = levsol_eig( fFit, vecGammaSeed, matHFit, matB, bMax, levPrm );
	%if ( stepSizeCoeff == 1.0 )
	%	vecZ = vecZFull;
	%elseif ( stepSizeCoeff < 1.0 )
	%	bMax = stepSizeCoeff*norm(vecZFull);
	%	vecZ = levsol_eig( fFit, vecGammaSeed, matHFit, matB, bMax, levPrm );
	%else
	%	stepSizeCoeff
	%	error( "stepSizeCoeff is inalid." );
	%endif
	%
	% 2022-12-29-2054: Sensible code.
	vecCap = max(abs(matVTDForFit),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	matB = diag(1.0./vecCap);
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	%msg( __FILE__, __LINE__, sprintf( "stepSizeCoeff = %0.3e", stepSizeCoeff ) );
	%msg( __FILE__, __LINE__, sprintf( "bMax = %0.3e", bMax ) );
	levPrm = [];
	assert( fFit >= 0.0 );
	vecZ = levsol_eig( fFit, vecGammaSeed, matHFit, matB, bMax, levPrm );
	if ( mygetfield( prm, "debugMode", false ) )
		vecDelta = matV*vecZ;
		vecCapFS = max(abs(matDForFit),[],2);
		vecCapFS += sqrt(eps)*max(vecCapFS);
		if ( abs(vecZ)./vecCap > bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
			[ vecDelta, vecCapFS, abs(vecDelta)./vecCapFS ]
			[ vecZ, vecCap, abs(vecZ)./vecCap ]
			bMax
		endif
		assert( abs(vecZ)./vecCap <= bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
	endif
	%
	% Note that fNew and vecGammaNew are merely estimates.
	vecYNew = vecYSeed + vecZ;
	fNew = fFit + vecYNew'*vecGammaFit + (vecYNew'*matHFit*vecYNew)/2.0;
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	normGammaNew = norm(vecGammaNew);
	%
	% We need to decide how to handle the changes to X and momentum that are outside of our subspace.
	% A few reasonable requirements...
	%  1. In the limit of the TR going to zero, the "jump" should leave us at the "seed".
	%  2. In the limit of fNew being zero, we should hit that point exactly with zero momentum.
	%  3. In the limit of taking the Newton step, the part of the gradient within the subspace should be zero.
	if ( fNew <= fSeed )
		alphaF = fNew / fSeed;
	else
		msg( __FILE__, __LINE__, "WARNING: fNew >= fSeed; this should never happen." );
		alphaF = 1.0;
	endif
	if ( normGammaNew < normGammaSeed )
		alphaG = normGammaNew / normGammaSeed;
	else
		alphaG = 1.0;
	endif
	% 2022-12-29-1740: These alpha are just reasonable guesses with some light testing.
	alphaXPerp = alphaG;
	alphaPPerp = alphaF;
	alphaGammaPerp = alphaG;
	%
	vecXNew = vecXAnchor + matV*vecYNew + vecXPerp*alphaXPerp;
	vecPNew = matV * ( coeffPGamma*vecGammaNew + vecGammaPerp*alphaGammaPerp ) + vecPPerp*alphaPPerp;
	%vecPNew = 0*vecXNew; % Is this the way to go???
return;
endfunction;


function [ vecXNew, vecPNew, jumpDat ] = __jump_simple_ineffective( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, trSize, prm )
	% DRaburn 2023-01-03:
	%  Because of issues with hessfit, decided to try a simpler fit.
	%  The fit itself works well, but this jump code doesn't.
	%  Not sure why, but, abandoning.
	%msg( __FILE__, __LINE__, "This jump does not work well when there is noise." );
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	%
	matDPreDrop = [ matX(:,1:indexAnchor-1), matX(:,indexAnchor+1:end) ] - vecXAnchor;
	matGPreDrop = [ matG(:,1:indexAnchor-1), matG(:,indexAnchor+1:end) ];
	rvecFPreDrop = [ rvecF(1:indexAnchor-1), rvecF(indexAnchor+1:end) ];
	orthoDropThresh = mygetfield( prm, "orthoDropThresh", 0.01 );
	[ matV, rvecDrop ] = utorthdrop( matDPreDrop, orthoDropThresh );
	%
	matDSans = matDPreDrop(:,~rvecDrop);
	matGSans = matGPreDrop(:,~rvecDrop);
	rvecFSans = rvecFPreDrop(~rvecDrop);
	jumpDat.sizeK = size(matDSans,2);
	% Now, we have our "sans (anchor)" and "anchor data".
	%
	matA = (matDSans'*matGSans) - (matDSans'*vecGAnchor); % Autobroadcast.
	matHFit = (matA'+matA)/2.0; % Simple. We could instead use "diagonal domination".
	fFit = fAnchor;
	%%%vecGammaFit = rvecFSans' - fAnchor - diag(matHFit)/2.0;
	vecGammaFit = matDSans'*vecGAnchor;
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matDSans'*vecGTrue;
		matHFit
		matVTHFSV = matDSans'*prm.matHSecret*matDSans
		matRes = matHFit - matVTHFSV
		rdH = reldiff( matHFit, matVTHFSV )
		vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
		rdG = reldiff( vecGammaFit, vecGammaTrue )
		fTrue
		fFit
		rdF = reldiff( fFit, fTrue )
	endif
	
	%
	%
	% We now have our (non-orthogonal) model...
	%   vecX = vecXAnchor +  matDeltaSans * vecY
	%   fModel = f0Fit + vecY'*vecGamma0Fit + (vecY'*matHFit*vecY)/2.0.
	vecD = sum(matDSans.^2, 1);
	%vecB = 1.0./( abs(vecD)+sqrt(eps)*max(abs(vecD)) );
	matB = diag(sqrt(vecD));
	%%%matB = symm(mtm(matDSans)^0.5);
	bMax = trSize;
	%bMax = [];
	levPrm = [];
	%
	vecXLaunch = vecXAnchor; % Projection of seed in to subspace may be better.
	fLaunch = fFit;
	vecGammaLaunch = vecGammaFit;
	vecY = levsol_eig( fLaunch, vecGammaLaunch, matHFit, matB, bMax, levPrm );
	%
	vecXNew = vecXLaunch + matDSans * vecY;
	vecPNew = vecPSeed - matV * ( matV' * vecPSeed ); % ... Or, something else?
	
	return;
	error( "END OF VALID CODE." );
	vecGammaNew = vecGammaLaunch + matHFit*vecY;
	%alphaG = norm(vecGammaNew)/norm(vecGammaLaunch);
	%vecPNew = vecPSeed*alphaG;
	%
	alphaG = norm(vecGammaNew)/norm(vecGammaLaunch);
	vecXPerp = NOTHING;
	%
	vecXNew = vecXLaunch + matDSans * vecY + vecXPerp*0
	[ norm(vecXNew-vecXLaunch), norm(matB*vecY), trSize ]
	vecPNew *= norm(vecGammaNew)/norm(vecGammaLaunch);
	%vecPNew = 0*vecXNew; % Is this the way to go???
return;
endfunction;


function [ vecXNew, vecPNew, jumpDat ] = __jump_simpleFitCts( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	
	% Fit is changed from original __jump_basicCts().
	%
	% Generate subspace basis matrix and drop columns.
	matD = matX - vecXAnchor; % Anchor column will necessarily be dropped.
	rvecDSq = sum( matD.^2, 1 );
	rvecDropBCMagnitude = ( rvecDSq < sqrt(eps)*max(rvecDSq) ); % Does almost nothing?!
	matDIntermed = matD(:,~rvecDropBCMagnitude);
	matGIntermed = matG(:,~rvecDropBCMagnitude);
	rvecFIntermed = rvecF(~rvecDropBCMagnitude);
	%
	basisDropThresh = mygetfield( prm, "basisDropThresh", 0.01 );
	[ matV, rvecDrop ] = utorthdrop( matDIntermed, basisDropThresh );
	jumpDat.sizeK = size(matV,2);
	if ( size(matV,2) == 0 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	matDSans = matDIntermed(:,~rvecDrop);
	matGSans = matGIntermed(:,~rvecDrop);
	rvecFSans = rvecFIntermed(~rvecDrop);
	%
	% y = V'*(x-x_anchor)
	% gamma = V'*g
	matY = matV'*matDSans;
	vecGammaAnchor = matV'*vecGAnchor;
	matGamma = matV'*matGSans;
	%
	% Generate fit.
	useTriuForFit = true;
	if (~useTriuForFit)
	useEigForFit = true;
	if (~useEigForFit)
	if ( rcond(matY*(matY')) < 100.0*eps )
		%matY
		size(matY)
		rc_mtm_yt = rcond(mtm(matY'))
		rc_mtm_y = rcond(mtm(matY))
		rd_triu = reldiff( matY, triu(matY))
		rc_mtm_triu_y = rcond(mtm(triu(matY)))
		rc_mtm_tril_yt = rcond(mtm(tril(matY')))
		rc_diag_y = rcond(diag(diag(matY)))
		rc_diag_mtm_y = rcond(diag(diag(mtm(matY))))
		rc_diag_mtm_yt = rcond(diag(diag(mtm(matY'))))
		%
		matYPrevIsh = matY(1:end-1,1:end-1);
		rc_mtm_ypt = rcond(mtm(matYPrevIsh'))
		rc_mtm_yp = rcond(mtm(matYPrevIsh))
	endif
	assert( rcond(matY'*matY) > 100.0*eps )
	assert( rcond(matY*(matY')) > 100.0*eps )
	matA = ( matGamma - vecGammaAnchor ) * (matY') / ( matY*(matY') ); % Autobroadcast.
	else
	[ matPsi, matLambda ] = eig(matY*(matY'));
	matA =  ( matGamma - vecGammaAnchor ) * (matY') * matPsi * inv(matLambda) * (matPsi');
	endif
	else
	matY_triu = triu(matY);
	matA = (matY_triu') \ (( matGamma - vecGammaAnchor)');
	endif
	matHFit = (matA'+matA)/2.0; % Alternatives are possible.
	fFit = fAnchor;
	vecGammaFit = vecGammaAnchor;
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matV'*vecGTrue;
		matHFit
		matVTHFSV = matV'*prm.matHSecret*matV
		matRes = matHFit - matVTHFSV
		rdH = reldiff( matHFit, matVTHFSV )
		vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
		rdG = reldiff( vecGammaFit, vecGammaTrue )
		fTrue
		fFit
		rdF = reldiff( fFit, fTrue )
	endif
	
	
	% Decompose (pre-jump) "seed" x.
	vecDSeed = vecXSeed - vecXAnchor;
	vecYSeed = matV'*vecDSeed;
	vecXPerp = vecDSeed - matV*vecYSeed;
	assert( reldiff( vecXSeed, vecXAnchor + (matV*vecYSeed) + vecXPerp ) <= sqrt(eps) );
	assert( norm(matV'*vecXPerp) <= sqrt(eps)*(norm(vecXAnchor)+norm(vecXSeed)) );
	%
	% These "seed" values are just estimates at the (subspace-projected) seed.
	fSeed = fFit + vecYSeed'*vecGammaFit + (vecYSeed'*matHFit*vecYSeed)/2.0;
	vecGammaSeed = vecGammaFit + matHFit*vecYSeed;
	normGammaSeed = norm(vecGammaSeed);
	%
	% Decompose (pre-jump) "seed" momentum.
	% Note that we actually do know how the full gradient varies throughout the subspace,
	%  not just how "gamma" (the part of the gradient that's in the subspace) varies;
	%  appropriate handling would essentially allow us to incorporate "vecPPerp" in to
	%  "coeffPGamma" and "vecGammaPerp".
	% Meh.
	vecT = matV'*vecPSeed;
	vecPPerp = vecPSeed - (matV*vecT);
	if ( 0.0 < normGammaSeed )
		coeffPGamma = (vecGammaSeed'*vecT) / (normGammaSeed^2);
	else
		coeffPGamma = 0.0;
		msg( __FILE__, __LINE__, "WARNING: 0.0 >= normGammaSeed; this should never happen." );
	endif
	vecGammaPerp = vecT - coeffPGamma*vecGammaSeed;
	assert( reldiff( matV*( coeffPGamma*vecGammaSeed + vecGammaPerp ) + vecPPerp, vecPSeed ) <= sqrt(eps) );
	assert( norm(matV'*vecPPerp) <= sqrt(eps)*norm(vecPSeed) );
	assert( abs(vecGammaPerp'*vecGammaSeed) <= sqrt(eps)*norm(vecPSeed) );
	%
	% Find point on Levenberg curve subject to trust region.
	%
	vecCap = max(abs(matY),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	matB = diag(1.0./vecCap);
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	%msg( __FILE__, __LINE__, sprintf( "stepSizeCoeff = %0.3e", stepSizeCoeff ) );
	%msg( __FILE__, __LINE__, sprintf( "bMax = %0.3e", bMax ) );
	levPrm = [];
	assert( fFit >= 0.0 );
	vecZ = levsol_eig( fFit, vecGammaSeed, matHFit, matB, bMax, levPrm );
	if ( mygetfield( prm, "debugMode", false ) )
		vecDelta = matV*vecZ;
		vecCapFS = max(abs(matDForFit),[],2);
		vecCapFS += sqrt(eps)*max(vecCapFS);
		if ( abs(vecZ)./vecCap > bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
			[ vecDelta, vecCapFS, abs(vecDelta)./vecCapFS ]
			[ vecZ, vecCap, abs(vecZ)./vecCap ]
			bMax
		endif
		assert( abs(vecZ)./vecCap <= bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
	endif
	%
	% Note that fNew and vecGammaNew are merely estimates.
	vecYNew = vecYSeed + vecZ;
	fNew = fFit + vecYNew'*vecGammaFit + (vecYNew'*matHFit*vecYNew)/2.0;
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	normGammaNew = norm(vecGammaNew);
	%
	% We need to decide how to handle the changes to X and momentum that are outside of our subspace.
	% A few reasonable requirements...
	%  1. In the limit of the TR going to zero, the "jump" should leave us at the "seed".
	%  2. In the limit of fNew being zero, we should hit that point exactly with zero momentum.
	%  3. In the limit of taking the Newton step, the part of the gradient within the subspace should be zero.
	if ( fNew <= fSeed )
		alphaF = fNew / fSeed;
	else
		msg( __FILE__, __LINE__, "WARNING: fNew >= fSeed; this should never happen." );
		alphaF = 1.0;
	endif
	if ( normGammaNew < normGammaSeed )
		alphaG = normGammaNew / normGammaSeed;
	else
		alphaG = 1.0;
	endif
	% 2022-12-29-1740: These alpha are just reasonable guesses with some light testing.
	alphaXPerp = alphaG;
	alphaPPerp = alphaF;
	alphaGammaPerp = alphaG;
	%
	vecXNew = vecXAnchor + matV*vecYNew + vecXPerp*alphaXPerp;
	vecPNew = matV * ( coeffPGamma*vecGammaNew + vecGammaPerp*alphaGammaPerp ) + vecPPerp*alphaPPerp;
	%vecPNew = 0*vecXNew; % Is this the way to go???
return;
endfunction;



function [ vecXNew, vecPNew, jumpDat ] = __jump_simpleFitReduced( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	% DRaburn 2023-01-03:
	%  This is per __jump_simpleFitCts but with stuff removed,
	%   to try to figure out why __jump_simple_ineffective was ultimately ineffective,
	%   and perhaps lead to elimination of the orthonormalization calculation.
	%  I've now trimmed it down to the presumed essentials (modulo the fact that stepSizeCoeff is wonky --
	%   should probably have both a "matD-based" and "external/BT-based" TR )
	%   and this still works.
	%  Still using orthonorm here vs _ineffetive, but there are also differences in effective TR.
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	
	% Fit is changed from original __jump_basicCts().
	%
	% Generate subspace basis matrix and drop columns.
	matD = matX - vecXAnchor; % Anchor column will necessarily be dropped.
	%
	% SKIP - Drop BC magnitude of columns of matD.
	%
	basisDropThresh = mygetfield( prm, "basisDropThresh", 0.01 );
	[ matV, rvecDrop ] = utorthdrop( matD, basisDropThresh );
	jumpDat.sizeK = size(matV,2);
	if ( size(matV,2) == 0 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	matDSans = matD(:,~rvecDrop);
	matGSans = matG(:,~rvecDrop);
	rvecFSans = rvecF(~rvecDrop);
	%
	% y = V'*(x-x_anchor)
	% gamma = V'*g
	matY = matV'*matDSans;
	vecGammaAnchor = matV'*vecGAnchor;
	matGamma = matV'*matGSans;
	%
	% Generate fit.
	matY_triu = triu(matY); % eig() on matY*(matY') could perhaps provide slightly more accurate results?
	matA = (matY_triu') \ (( matGamma - vecGammaAnchor)');
	matHFit = (matA'+matA)/2.0; % Alternatives are possible.
	fFit = fAnchor;
	vecGammaFit = vecGammaAnchor;
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matV'*vecGTrue;
		matHFit
		matVTHFSV = matV'*prm.matHSecret*matV
		matRes = matHFit - matVTHFSV
		rdH = reldiff( matHFit, matVTHFSV )
		vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
		rdG = reldiff( vecGammaFit, vecGammaTrue )
		fTrue
		fFit
		rdF = reldiff( fFit, fTrue )
	endif
	
	%
	% SKIP - Find nearest point to vecXSeed to use as launch.
	%
	%
	% Find point on Levenberg curve subject to trust region.
	%
	vecYLaunch = zeros(size(vecGammaFit));
	fLaunch = fFit;
	vecGammaLaunch = vecGammaFit;
	vecCap = max(abs(matY),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	matB = diag(1.0./vecCap);
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	levPrm = [];
	vecZ = levsol_eig( fLaunch, vecGammaLaunch, matHFit, matB, bMax, levPrm );
	if ( mygetfield( prm, "debugMode", false ) )
		vecDelta = matV*vecZ;
		vecCapFS = max(abs(matDForFit),[],2);
		vecCapFS += sqrt(eps)*max(vecCapFS);
		if ( abs(vecZ)./vecCap > bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
			[ vecDelta, vecCapFS, abs(vecDelta)./vecCapFS ]
			[ vecZ, vecCap, abs(vecZ)./vecCap ]
			bMax
		endif
		assert( abs(vecZ)./vecCap <= bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
	endif
	%
	% Note that fNew and vecGammaNew are merely estimates.
	vecYNew = vecYLaunch + vecZ;
	fNew = fFit + vecYNew'*vecGammaFit + (vecYNew'*matHFit*vecYNew)/2.0;
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	%
	% SKIP - Decomposition of vecXSeed and vecPSeed and analysis
	%  to get improved vecXNew and vecPNew, as well as make output "new" values
	%  approach input "seed" values in the limit of step size -> 0.
	vecXNew = vecXAnchor + matV*vecYNew;
	vecPNew = vecPSeed * norm(vecGammaNew) / norm(vecGammaFit);
return;
endfunction;


function [ vecXNew, vecPNew, jumpDat ] = __jump_simple_redux( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	% DRaburn 2023-01-03:
	%  Another attempt at not using orthogonalization.
	%  Does not yet work; try matching TR in __jump_simpleFitReduced exactly?
	%  Maybe fit is less accurate because of non-ortho subspace???
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	%
	matD = matX - vecXAnchor;
	orthoDropThresh = mygetfield( prm, "orthoDropThresh", 0.01 );
	[ foo, rvecDrop ] = utorthdrop( matD, orthoDropThresh );
	matDSans = matD(:,~rvecDrop);
	matGSans = matG(:,~rvecDrop);
	rvecFSans = rvecF(~rvecDrop);
	jumpDat.sizeK = size(matDSans,2);
	%
	matA = (matDSans'*matGSans) - (matDSans'*vecGAnchor); % Autobroadcast.
	matHFit = (matA'+matA)/2.0; % Simple. We could instead use "diagonal domination".
	fFit = fAnchor;
	vecGammaFit = matDSans'*vecGAnchor;
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matDSans'*vecGTrue;
		matHFit
		matVTHFSV = matDSans'*prm.matHSecret*matDSans
		matRes = matHFit - matVTHFSV
		rdH = reldiff( matHFit, matVTHFSV )
		vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
		rdG = reldiff( vecGammaFit, vecGammaTrue )
		fTrue
		fFit
		rdF = reldiff( fFit, fTrue )
	endif
	
	%
	matB = [];
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	levPrm = [];
	%
	vecYLaunch = zeros(size(vecGammaFit));
	fLaunch = fFit;
	vecGammaLaunch = vecGammaFit;
	vecZ = levsol_eig( fLaunch, vecGammaLaunch, matHFit, matB, bMax, levPrm );
	vecYNew = vecYLaunch + vecZ;
	%
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	%
	vecXNew = vecXAnchor + matDSans * vecYNew;
	vecPNew = vecPSeed * norm( vecGammaNew ) / norm( vecGammaLaunch );
return;
endfunction;


function [ vecXNew, vecPNew, jumpDat ] = __jump_january( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	% DRaburn 2023-01-04:
	%  About ready to be done with this solver.
	%  But, let's try out a few things first.
	%  Copied from __jump_simpleFitReduced.
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	
	% Fit is changed from original __jump_basicCts().
	%
	% Generate subspace basis matrix and drop columns.
	matD = matX - vecXAnchor; % Anchor column will necessarily be dropped.
	%
	% SKIP - Drop BC magnitude of columns of matD.
	%
	basisDropThresh = mygetfield( prm, "basisDropThresh", 0.01 );
	[ matV, rvecDrop ] = utorthdrop( matD, basisDropThresh );
	jumpDat.sizeK = size(matV,2);
	if ( size(matV,2) == 0 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	matDSans = matD(:,~rvecDrop);
	matGSans = matG(:,~rvecDrop);
	rvecFSans = rvecF(~rvecDrop);
	%
	% y = V'*(x-x_anchor)
	% gamma = V'*g
	matY = triu(matV'*matDSans);
	vecGammaAnchor = matV'*vecGAnchor;
	matGamma = matV'*matGSans;
	%
	% Generate fit.
	fFit = fAnchor;
	vecGammaFit = vecGammaAnchor;
	matA = (matY') \ (( matGamma - vecGammaAnchor)');
	matHFit = (matA'+matA)/2.0;
	% DRaburn 2023-1-04:
	% Alternatives approaches to matHFit are possible,
	%  as, considering rvecF, we have two measures of every element of matHFit.
	% However, my attempts did not work particularly well and added to the time.
	% See hessfit_simple_posdefy.m, for example.
	
	doComparison = false;
	if (doComparison)
		[ fTrue, vecGTrue ] = prm.funchFG_noiseless( vecXAnchor );
		vecGammaTrue = matV'*vecGTrue;
		matVTHFSV = matV'*prm.matHSecret*matV;
		rdH = reldiff( matHFit, matVTHFSV );
		rdG = reldiff( vecGammaFit, vecGammaTrue );
		rdF = reldiff( fFit, fTrue );
		%
		dumpInfo = true;
		if (dumpInfo)
			matHFit
			matVTHFSV
			matRes = matHFit - matVTHFSV
			vecGammaCompare = [ vecGammaFit, vecGammaTrue, vecGammaFit - vecGammaTrue ]
			fCompare = [ fFit, fTrue, fFit - fTrue ]
			rdVals = [ rdF, rdG, rdH ]
		endif
		if ( rdH > 0.1 )
			error( "*** THE BAD THING HAPPENED ***" );
		endif
	endif
	
	%
	% SKIP - Find nearest point to vecXSeed to use as launch.
	%
	%
	% Find point on Levenberg curve subject to trust region.
	%
	vecYLaunch = zeros(size(vecGammaFit));
	fLaunch = fFit;
	vecGammaLaunch = vecGammaFit;
	vecCap = max(abs(matY),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	matB = diag(1.0./vecCap);
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	levPrm = [];
	vecZ = levsol_eig( fLaunch, vecGammaLaunch, matHFit, matB, bMax, levPrm );
	if ( mygetfield( prm, "debugMode", false ) )
		vecDelta = matV*vecZ;
		vecCapFS = max(abs(matDForFit),[],2);
		vecCapFS += sqrt(eps)*max(vecCapFS);
		if ( abs(vecZ)./vecCap > bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
			[ vecDelta, vecCapFS, abs(vecDelta)./vecCapFS ]
			[ vecZ, vecCap, abs(vecZ)./vecCap ]
			bMax
		endif
		assert( abs(vecZ)./vecCap <= bMax*1.1 || abs(vecDelta)./vecCapFS > bMax*1.1 )
	endif
	%
	% Note that fNew and vecGammaNew are merely estimates.
	vecYNew = vecYLaunch + vecZ;
	fNew = fFit + vecYNew'*vecGammaFit + (vecYNew'*matHFit*vecYNew)/2.0;
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	%
	% SKIP - Decomposition of vecXSeed and vecPSeed and analysis
	%  to get improved vecXNew and vecPNew, as well as make output "new" values
	%  approach input "seed" values in the limit of step size -> 0.
	vecXNew = vecXAnchor + matV*vecYNew;
	vecPNew = vecPSeed * norm(vecGammaNew) / norm(vecGammaFit);
return;
endfunction;


function [ vecXNew, vecPNew, jumpDat ] = __jump_januarySeedy( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
	% DRaburn 2023-01-04:
	%  Copied from __jump_january().
	%  Let's add "seed" handling per __jump_basicCts():
	%   consideration of subspace and being continuous.
	%  ( Side note: basicCts no longer seems very basic. )
	vecXNew = [];
	vecPNew = [];
	jumpDat = [];
	jumpDat.sizeK = 0; % Unless overwritten.
	jumpDat.rvecUseForFit = []; % Unless overwritten.
	%
	if ( size(matX,2) <= 1 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	%
	% Select anchor for subspace.
	[ fAnchor, indexAnchor ] = min(rvecF);
	vecXAnchor = matX(:,indexAnchor);
	vecGAnchor = matG(:,indexAnchor);
	%
	% Generate subspace basis matrix and drop columns.
	matD = matX - vecXAnchor; % Anchor column will necessarily be dropped.
	% We could consider dropping columns because of the magnitude of the corresponding vector of matD is too small.
	% This is done in <some earlier __jump version>.
	basisDropThresh = mygetfield( prm, "basisDropThresh", 0.01 );
	[ matV, rvecDrop ] = utorthdrop( matD, basisDropThresh );
	jumpDat.sizeK = size(matV,2);
	if ( size(matV,2) == 0 )
		vecXNew = vecXSeed;
		vecPNew = vecPSeed;
		return;
	endif
	matDSans = matD(:,~rvecDrop);
	matGSans = matG(:,~rvecDrop);
	rvecFSans = rvecF(~rvecDrop);
	%
	% y = V'*(x-x_anchor)
	% gamma = V'*g
	matY = triu(matV'*matDSans);
	vecGammaAnchor = matV'*vecGAnchor;
	matGamma = matV'*matGSans;
	%
	% Generate fit.
	fFit = fAnchor;
	vecGammaFit = vecGammaAnchor;
	matA = (matY') \ (( matGamma - vecGammaAnchor)');
	matHFit = (matA'+matA)/2.0;
	% DRaburn 2023-01-04:
	% Alternatives approaches to matHFit are possible. See other versions of jump code.
	% Omitting fit comparision.
	%
	
	% Decompose (pre-jump) "seed" vecX.
	vecDSeed = vecXSeed - vecXAnchor;
	vecYSeed = matV'*vecDSeed;
	vecXPerp = vecDSeed - matV*vecYSeed;
	assert( reldiff( vecXSeed, vecXAnchor + (matV*vecYSeed) + vecXPerp ) <= sqrt(eps) );
	assert( norm(matV'*vecXPerp) <= sqrt(eps)*(norm(vecXAnchor)+norm(vecXSeed)) );
	%
	% Note: The following "seed" values are *estimates* at the (subspace-projected) seed.
	fSeed = fFit + vecYSeed'*vecGammaFit + (vecYSeed'*matHFit*vecYSeed)/2.0;
	vecGammaSeed = vecGammaFit + matHFit*vecYSeed;
	normGammaSeed = norm(vecGammaSeed);
	if ( normGammaSeed <= 0.0 )
		msg( __FILE__, __LINE__, "WARNING: norm(vecGammaSeed) <= 0.0; this is unexpected!!!" );
	endif
	%
	% Decompose (pre-jump) "seed" momentum.
	% Note that we actually *do* have information about how the full gradient varies throughout the subspace,
	%  not just how "gamma" (the part of the gradient that's in the subspace) varies;
	%  appropriate handling would essentially allow us to incorporate "vecPPerp" into "coeffPGamma" and "vecGammaPerp".
	vecT = matV'*vecPSeed;
	vecPPerp = vecPSeed - (matV*vecT);
	assert( reldiff( vecPSeed, matV*vecT + vecPPerp ) <= sqrt(eps) );
	% Further, we decompose the temporary vector...
	if ( normGammaSeed > 0.0 )
		coeffPGamma = (vecGammaSeed'*vecT) / (normGammaSeed^2);
	else
		coeffPGamma = 0.0;
	endif
	vecGammaPerp = vecT - coeffPGamma*vecGammaSeed;
	assert( reldiff( vecPSeed, matV*( coeffPGamma*vecGammaSeed + vecGammaPerp ) + vecPPerp ) <= sqrt(eps) );
	assert( norm(matV'*vecPPerp) <= sqrt(eps)*norm(vecPSeed) );
	assert( abs(vecGammaPerp'*vecGammaSeed) <= sqrt(eps)*norm(vecPSeed) );
	
	%
	% Set up Levenberg curve.
	vecYLaunch = vecYSeed;
	fLaunch = fFit + vecGammaFit'*vecYLaunch + (vecYLaunch'*matHFit*vecYLaunch)/2.0;
	vecGammaLaunch = vecGammaFit + matHFit*vecYLaunch;
	%
	vecCap = max(abs(matY),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	matB = diag(1.0./vecCap);
	trCoeff = mygetfield( prm, "trCoeff", 3.0 );
	bMax = trCoeff * stepSizeCoeff;
	levPrm = [];
	%
	vecZ = levsol_eig( fLaunch, vecGammaLaunch, matHFit, matB, bMax, levPrm );
	% Omitting checks on "cap/matB".
	%
	% Note that fNew and vecGammaNew are merely estimates.
	vecYNew = vecYLaunch + vecZ;
	fNew = fFit + vecYNew'*vecGammaFit + (vecYNew'*matHFit*vecYNew)/2.0;
	vecGammaNew = vecGammaFit + matHFit*vecYNew;
	
	%
	% We need to decide how to handle the changes to X and momentum that are outside of our subspace.
	% A few reasonable requirements...
	%  1. In the limit of the TR going to zero, the "jump" should leave us at the "seed".
	%  2. In the limit of fNew being zero, we should hit that point exactly with zero momentum.
	%  3. In the limit of taking the Newton step, the part of the gradient within the subspace should be zero.
	normGammaNew = norm(vecGammaNew);
	if ( fNew <= fSeed )
		alphaF = fNew / fSeed;
	else
		msg( __FILE__, __LINE__, "WARNING: fNew >= fSeed; this should never happen." );
		alphaF = 1.0;
	endif
	if ( normGammaNew < normGammaSeed )
		alphaG = normGammaNew / normGammaSeed;
	else
		alphaG = 1.0;
	endif
	% 2022-12-29-1740: These alpha are just reasonable guesses with some light testing.
	alphaXPerp = alphaG;
	alphaPPerp = alphaF;
	alphaGammaPerp = alphaG;
	%
	vecXNew = vecXAnchor + matV*vecYNew + vecXPerp*alphaXPerp;
	vecPNew = matV * ( coeffPGamma*vecGammaNew + vecGammaPerp*alphaGammaPerp ) + vecPPerp*alphaPPerp;
	%vecPNew = 0*vecXNew; % Is this the way to go???
	
return;
endfunction;
