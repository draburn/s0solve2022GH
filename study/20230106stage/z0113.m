clear;
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
	return;
endif
stopsig_lastTime = time();
%[ funchFG, vecX0, prm, initDat ] = z__init0113();
%[ funchFG, vecX0, prm, initDat ] = z__init0125danmoqnjpy();
[ funchFG, vecX0, prm, initDat ] = z__init20250501();
sizeX = size(vecX0,1);
%
vecX = vecX0;
vecP = zeros(sizeX,1);
vecXSeed = vecX;
vecPSeed = vecP;
%
best_vecX = vecX0;
best_vecG = [];
best_f = [];
best_fVar = [];
best_vecXHarvest = [];
best_vecPHarvest = [];
minf_vecX = vecX0;
minf_vecG = [];
minf_f = [];
minf_fVar = [];
%
startTime = time();
iterCount = 0;
fevalCount = 0;
badCount = 0;
%
proglog_lastTime = time();
proglog_stopTime = time();
%
running_fevalCount = 0;
running_fTot = 0.0;
running_xtgTot = 0.0;
running_vecGTot = zeros(sizeX,1);
running_vecXTot = zeros(sizeX,1);
running_fSqTot = 0.0;
running_xtgSqTot = 0.0;
running_vecGSqTot = zeros(sizeX,1);
running_vecXSqTot = zeros(sizeX,1);
%
% Only really used by QNJ
record_matX = [];
record_matG = [];
record_rvecF = [];
record_rvecW = [];
record_indexToDrop = []; % Mis-named?
qnj_matV = [];
qnj_sMax = prm.qnj_sMaxInit;
qnj_dMax = prm.qnj_dMaxInit;
qnj_s = [];
qnj_vecDelta = [];
%
while ( 1 )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CONVENTIONAL GRADIENT DESCENT + STUFF
	%
	[ f, vecG ] = funchFG( vecX );
	fevalCount++;
	%
	xtg = vecX'*vecG;
	running_fevalCount++;
	running_fTot += f;
	running_xtgTot += xtg;
	running_vecGTot += vecG;
	running_vecXTot += vecX;
	running_fSqTot += f^2;
	running_xtgSqTot += xtg^2;
	running_vecGSqTot += vecG.^2;
	running_vecXSqTot += vecX.^2;
	%
	vecP = ( prm.momentumFactor * vecP ) - ( prm.learningRate * vecG );
	vecX += vecP;
	%
	if ( f > prm.fBail )
		msg( __FILE__, __LINE__, "IMPOSED STOP: f > prm.fBail. This strongly indicates divergence." );
		break;
	elseif ( prm.fevalLimit > 0 && fevalCount >= prm.fevalLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: fevalCount >= prm.fevalLimit." );
		break;
	elseif ( prm.timeLimit > 0 && time() - startTime >= prm.timeLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: time() - startTime >= prm.timeLimit." );
		break;
	elseif ( prm.stopSignalCheckInterval >= 0.0 && time() - stopsig_lastTime >= prm.stopSignalCheckInterval )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		stopsig_lastTime = time();
	endif
	%
	if ( running_fevalCount < prm.numFevalPerSuperPt )
		continue
	endif
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SUPER-POINT ANALYSIS
	%
	vecXHarvest = vecX;
	vecPHarvest = vecP;
	iterCount++;
	%
	superPt_vecX = running_vecXTot / running_fevalCount;
	superPt_vecG = running_vecGTot / running_fevalCount;
	superPt_xtgAvg = running_xtgTot / running_fevalCount;
	superPt_fAvg = running_fTot / running_fevalCount;
	superPt_f = superPt_fAvg - (( superPt_xtgAvg - (superPt_vecX'*superPt_vecG) )/2.0);
	superPt_w = running_fevalCount;
	superPt_vecXSqVar = (running_vecXSqTot / running_fevalCount) - (superPt_vecX.^2);
	superPt_vecGSqVar = (running_vecGSqTot / running_fevalCount) - (superPt_vecG.^2);
	superPt_xtgSqVar = (running_xtgSqTot / running_fevalCount) - (superPt_xtgAvg^2);
	superPt_fSqVar = (running_fSqTot / running_fevalCount) - (superPt_fAvg^2);
	superPt_fVar = real(sqrt(superPt_fSqVar)); % Note part of this variation is due to variation of x.
	superPt_gVar = real(sqrt(sum(superPt_vecGSqVar))); % Note part of this variation is due to variation of x.
	%
	running_fevalCount = 0;
	running_fTot = 0.0;
	running_xtgTot = 0.0;
	running_vecGTot = zeros(sizeX,1);
	running_vecXTot = zeros(sizeX,1);
	running_fSqTot = 0.0;
	running_xtgSqTot = 0.0;
	running_vecGSqTot = zeros(sizeX,1);
	running_vecXSqTot = zeros(sizeX,1);
	%
	newIsBest = false; % Unless...
	newIsMinf = false; % Unless...
	if ( isempty(minf_f) )
		newIsBest = true;
		newIsMinf = true;
	elseif ( superPt_f < minf_f  )
		newIsBest = true;
		newIsMinf = true;
	elseif ( (superPt_f <= minf_f + prm.bestFVarCoeffA*superPt_fVar + prm.bestFVarCoeffB*minf_fVar) ...
	  && (norm(superPt_vecG) < norm(best_vecG)) )
		newIsBest = true;
	endif
	if ( newIsBest )
		best_vecX = superPt_vecX;
		best_vecG = superPt_vecG;
		best_f = superPt_f;
		best_fVar = superPt_fVar;
		best_vecXHarvest = vecXHarvest;
		best_vecPHarvest = vecPHarvest;
	else
		badCount++;
	endif
	if ( newIsMinf )
		minf_vecX = superPt_vecX;
		minf_vecG = superPt_vecG;
		minf_f = superPt_f;
		minf_fVar = superPt_fVar;
	endif
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CHECK STOP CRIT
	%
	if ( norm( superPt_vecG ) <= prm.gTol )
		msg( __FILE__, __LINE__, "SUCCESS: norm( superPt_vecG ) <= prm.gTol." );
		break;
	elseif ( superPt_f <= prm.fTol )
		msg( __FILE__, __LINE__, "SUCCESS: superPt_f <= prm.fTol." );
		break;
	%elseif ( superPt_f < superPt_fPrev && superPt_fPrev - superPt_f <= prm.fTol )
	%	msg( __FILE__, __LINE__, "IMPOSED STOP: superPt_f < superPt_fPrev && superPt_fPrev - superPt_f <= prm.fTol." );
	%	break;
	elseif ( norm( vecXHarvest - vecXSeed ) <= prm.xTol )
		msg( __FILE__, __LINE__, "IMPOSED STOP: norm( vecXHarvest - vecXSeed ) <= prm.xTol." );
		break;
	elseif ( prm.iterLimit >= 0 && iterCount >= prm.iterLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterLimit." );
		break;
	endif
	if ( time() > proglog_lastTime + prm.progressReportInterval )
		z__proglog0113;
	 	proglog_lastTime = time();
	endif
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PREPARE FOR NEXT ITER (ESP. NON-QNJ)
	%
	if ( size(record_matX,2) < prm.maxNumRecords )
		record_matX = [ superPt_vecX, record_matX ];
		record_matG = [ superPt_vecG, record_matG ];
		record_rvecF = [ superPt_f, record_rvecF ];
		record_rvecW = [ superPt_w, record_rvecW ];
	elseif ( isempty(record_indexToDrop) )
		record_matX = [ superPt_vecX, record_matX(:,1:end-1) ];
		record_matG = [ superPt_vecG, record_matG(:,1:end-1) ];
		record_rvecF = [ superPt_f, record_rvecF(1:end-1) ];
		record_rvecW = [ superPt_w, record_rvecW(1:end-1) ];
	else
		matX = [ superPt_vecX, matX(:,1:record_indexToDrop-1,record_indexToDrop+1:end) ];
		matG = [ superPt_vecG, matG(:,1:record_indexToDrop-1,record_indexToDrop+1:end) ];
		rvecF = [ superPt_f, rvecF(1:record_indexToDrop-1,record_indexToDrop+1:end) ];
		rvecW = [ superPt_w, rvecW(1:record_indexToDrop-1,record_indexToDrop+1:end) ];
	endif
	%
	% These results are likely to be modified by QNJ, if used.
	vecX = vecXHarvest; % what actually gets used as the seed.
	vecP = vecPHarvest;
	vecXSeed = vecX; % store, for reference.
	vecPSeed = vecP;
	qnj_matV = [];
	
	doEarlyBasisGeneration = true;
	if (doEarlyBasisGeneration)
		% Doing this to provide comparision to python code basis generation.
		vecXAnchor = best_vecX;
		vecGAnchor = best_vecG;
		fAnchor = best_f;
		matD = record_matX - vecXAnchor;
		%record_matG
		[ qnj_matV, rvecBasisDrop ] = utorthdrop( matD, prm.qnj_basisDropThresh );
		%qnj_matV
		if ( size(qnj_matV,2) < 1 )
			continue;
		endif
		matY = triu( qnj_matV' * matD(:,~rvecBasisDrop) );
		%msg( __FILE__, __LINE__, "matY = ..." ); matY
		matGamma = qnj_matV' * record_matG(:,~rvecBasisDrop);
		vecGammaAnchor = qnj_matV' * vecGAnchor;
	endif
	
	%
	if ( ~prm.useQNJ )
		continue;
	endif
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO QNJ WORK
	%
	% Or, maybe not?
	if ( prm.qnj_forceGradAfterBad && ~newIsBest && ~isempty(best_vecXHarvest) && ~isempty(best_vecPHarvest) )
		vecX = best_vecXHarvest;
		vecP = best_vecPHarvest;
		vecXSeed = vecX; % store, for reference.
		vecPSeed = vecP;
		continue;
		% I suppose this means that, if the grad steps continue to be bad, we'll never attempt QNJ. Fine.
	endif
	
	if (~doEarlyBasisGeneration)
	% Generate basis.
	% 2023-01-13: I know of no better way to do this, quality-wise.
	vecXAnchor = best_vecX;
	vecGAnchor = best_vecG;
	fAnchor = best_f;
	matD = record_matX - vecXAnchor;
	[ qnj_matV, rvecBasisDrop ] = utorthdrop( matD, prm.qnj_basisDropThresh );
	if ( size(qnj_matV,2) < 1 )
		continue;
	endif
	matY = triu( qnj_matV' * matD(:,~rvecBasisDrop) );
	matGamma = qnj_matV' * record_matG(:,~rvecBasisDrop);
	vecGammaAnchor = qnj_matV' * vecGAnchor;
	endif
	
	%
	% Generate fit.
	% 2023-01-13: This is a bit simplistic but reasonable.
	fFit = fAnchor;
	vecGammaFit = vecGammaAnchor;
	matA = (matY') \ (( matGamma - vecGammaAnchor )');
	matHFit = (matA'+matA)/2.0;
	%
	% Update trust region.
	% 2023-01-13: This deserves testing.
	if ( ~isempty(qnj_vecDelta) )
		if ( newIsBest )
			% 2023-01-15: 'limitOnGood' was formerly called 'pre-limit'.
			if ( prm.qnj_sMax_limitOnGood )
				qnj_sMax = prm.qnj_sMaxFT * qnj_s;
			elseif ( ~isempty(qnj_sMax) )
				qnj_sMax = max([ prm.qnj_sMaxFT * qnj_s, qnj_sMax ]);
			endif
			if ( prm.qnj_dMax_limitOnGood )
				qnj_dMax = prm.qnj_dMaxFT * norm(qnj_vecDelta);
			elseif ( ~isempty(qnj_dMax) )
				qnj_dMax = max([ prm.qnj_dMaxFT * norm(qnj_vecDelta), qnj_dMax ]);
			endif
		else
			qnj_sMax = prm.qnj_sMaxBT * qnj_s;
			qnj_dMax = prm.qnj_dMaxBT * norm(qnj_vecDelta);
		endif
	endif
	qnj_sMax = mycap( qnj_sMax, prm.qnj_sMaxLo, prm.qnj_sMaxHi );
	qnj_dMax = mycap( qnj_dMax, prm.qnj_dMaxLo, prm.qnj_dMaxHi );
	vecCap = max(abs(matY),[],2);
	vecCap += sqrt(eps)*max(vecCap);
	assert( 0.0 < min(vecCap) );
	vecS = 1.0 ./ vecCap;
	% 2023-01-15: We could also consider a limit on fMin and fModMin.
	%
	%
	% Calculate the next guess.
	if (prm.qnj_useCtsFromHarvest)
		% Decompose vecXLaunch.
		vecXLaunch = best_vecXHarvest;
		vecDLaunch = vecXLaunch - vecXAnchor;
		vecYLaunch = qnj_matV' * vecDLaunch;
		vecXPerp = vecDLaunch - ( qnj_matV * vecYLaunch );
		assert( reldiff( vecXLaunch, vecXAnchor + (qnj_matV*vecYLaunch) + vecXPerp ) <= sqrt(eps) );
		assert( norm(qnj_matV'*vecXPerp) <= sqrt(eps)*(norm(vecXAnchor)+norm(vecXLaunch)) );
		%
		% Calculate fLaunch and vecGammaLaunch as close as possible to vecXLaunch.
		% We need matHMod.
		% Calling levsol is computationally wasteful, but POITROME.
		[ aux_vecZNewton, aux_stepDat ] = levsol0111( fFit, vecGammaFit, matHFit, vecS );
		fLaunch = fFit + (vecYLaunch'*vecGammaFit) + ((vecYLaunch'*aux_stepDat.matHMod*vecYLaunch)/2.0);
		vecGammaLaunch = vecGammaFit + (aux_stepDat.matHMod*vecYLaunch);
		%
		% Decompose vecPLaunch.
		vecPLaunch = best_vecPHarvest;
		vecT = qnj_matV' * vecPLaunch;
		vecPPerp = vecPLaunch - ( qnj_matV * vecT );
		assert( norm(vecGammaLaunch) > 0.0 );
		coeffPG = (vecGammaLaunch'*vecT)/sumsq(vecGammaLaunch);
		vecGammaPerp = vecT - ( coeffPG * vecGammaLaunch );
		assert( reldiff( vecPLaunch, (qnj_matV*( (coeffPG*vecGammaLaunch) + vecGammaPerp )) + vecPPerp ) <= sqrt(eps) );
		assert( norm(qnj_matV'*vecPPerp) <= sqrt(eps)*(norm(vecPLaunch)+norm(vecGammaPerp)+abs(coeffPG)*norm(vecGammaLaunch)) );
		if ( coeffPG > 0.0 )
			%msg( __FILE__, __LINE__, sprintf( "WARNING: Momentum was going uphill (%0.3E)", coeffPG ) );
			coeffPG = 0.0;
		endif
		%
		% Calculate step.
		stepPrm = [];
		[ vecZ, stepDat ] = levsol0111( fLaunch, vecGammaLaunch, matHFit, vecS, qnj_sMax, qnj_dMax, stepPrm );
		vecYNew = vecYLaunch + vecZ;
		fNew = fLaunch + (vecZ'*vecGammaLaunch) + ((vecZ'*stepDat.matHMod*vecZ)/2.0);
		vecGammaNew = vecGammaLaunch + (stepDat.matHMod*vecZ);
		assert( fLaunch > 0.0 );
		assert( fNew < fLaunch );
		%assert( fNew*(1.0-100.0*eps) <= fLaunch*(1.0+100.0*eps) );
		alphaF = fNew / fLaunch;
		assert( norm(vecS.\vecGammaNew) < norm(vecS.\vecGammaLaunch) );
		%assert( norm(vecS.\vecGammaNew)*(1.0-100.0*eps) < norm(vecS.\vecGammaLaunch)*(1.0+100.0*eps) );
		alphaG = norm(vecS.\vecGammaNew) / norm(vecS.\vecGammaLaunch);
		%
		vecX = vecXAnchor + (qnj_matV * vecYNew) + (alphaG * vecXPerp);
		vecP = (qnj_matV*( (coeffPG*vecGammaNew) + (alphaG*vecGammaPerp) )) + (alphaF*vecPPerp);
		%
		qnj_vecDelta = vecX - vecXHarvest;
		qnj_s = norm( vecS .* vecZ );
	else
	% 2023-01-13: This is placeholder; compare to study/20221216sxsolve/sgsolve.m.
	vecYLaunch = zeros(size(vecGammaFit));
	fLaunch = fFit + vecYLaunch'*vecGammaFit + (vecYLaunch'*matHFit*vecYLaunch)/2.0;
	vecGammaLaunch = vecGammaFit + ( matHFit * vecYLaunch );
	stepPrm = [];
	[ vecZ, stepDat ] = levsol0111( fLaunch, vecGammaLaunch, matHFit, vecS, qnj_sMax, qnj_dMax, stepPrm );
	%[ vecZ, stepDat ] = eigfloorsol0111( fLaunch, vecGammaLaunch, matHFit, vecS, qnj_sMax, qnj_dMax, stepPrm );
	assert( ~isempty(vecZ) );
	vecYNew = vecYLaunch + vecZ;
	vecGammaNew = stepDat.vecGModPred;
	gammaRatio = norm(vecS.\vecGammaNew) / norm(vecS.\vecGammaLaunch); % Need scaling to be monotonic.
	assert( gammaRatio <= 1.0 );
	%qnj_vecDelta = (qnj_matV * vecYNew) + vecXAnchor - vecX; HORRID
	%qnj_vecDelta = (qnj_matV * vecYNew) + (1.0-gammaRatio)*( vecXAnchor - vecX ); BAD
	qnj_vecDelta = (qnj_matV * vecYNew); % GOOD. Why????
	qnj_s = norm( vecS .* vecZ );
	vecX += qnj_vecDelta;
	vecP -= (1.0-gammaRatio)*vecP;
	endif
	%
	% Decide which record to drop... or merge???
	% 2023-01-13: This is crude.
	record_indexToDrop = [];
	%
	% And, we are outt'a here!
	vecXSeed = vecX;
	vecPSeed = vecP;
endwhile
vecXF = best_vecX;
z__proglog0113;
z__fin0113;
