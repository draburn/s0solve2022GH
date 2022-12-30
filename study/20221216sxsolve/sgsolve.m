function [ vecX, retCode, datOut ] = sgsolve( funchFG, init_vecX, prm=[] )
	msg( __FILE__, __LINE__, "Need:" );
	msg( __FILE__, __LINE__, "  ( Testing? Refactor? )" );
	msg( __FILE__, __LINE__, "Want Eventually:" );
	msg( __FILE__, __LINE__, "  * Reduce repeated manipulation of record data (matX, matG)." );
	msg( __FILE__, __LINE__, "Maybe Consider:" );
	msg( __FILE__, __LINE__, "  * Implement two- and three-pass hessfit, possibly using non-orthogonal basis." );
	msg( __FILE__, __LINE__, "  * Model variation of gradient outside V when jump." );
	msg( __FILE__, __LINE__, "  * Find analytic justification for 'alpha' values when jump, for seed data outside subspace." );
	msg( __FILE__, __LINE__, "  * Enable and test non-orthogonal basis matrix (matD except for drop) for hessfit; reduce FP issues?" );
	msg( __FILE__, __LINE__, "  * Intelligently curate record data (matX, matG), balancing num fevals and own work." );
	msg( __FILE__, __LINE__, "  * Further test TR." );
	msg( __FILE__, __LINE__, "Potential Optimizations:" );
	msg( __FILE__, __LINE__, "  * Improved 'momentum' via estimating gradient at 'sprout' point." );
	msg( __FILE__, __LINE__, "  * Super-point data from just start & end of 'super-point/meander' instead of per-feval." );
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
		stepSizeCoeff = 1.0;
		[ seed_vecX, seed_vecP, jumpDat ] = __jump_basicCts( seed_vecX, seed_vecP, matX, matG, rvecF, rvecW, stepSizeCoeff, prm );
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
	datOut.wSPt = wSum;
return;
endfunction


function [ vecXNew, vecPNew, jumpDat ] = __jump_basicCts( vecXSeed, vecPSeed, matX, matG, rvecF, rvecW, stepSizeCoeff, prm )
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
	[ fFit, vecGammaFit, matHFit ] = hessfit( matVTDForFit, rvecFForFit, matVTGForFit, fitPrm );
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
	levPrm = [];
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
return;
endfunction;
