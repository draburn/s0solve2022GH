clear;
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
	return;
endif
stopsig_lastTime = time();
[ funchFG, vecX0, prm, initDat ] = z__init0113();
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
minf_vecX = vecX0;
minf_vecG = [];
minf_f = [];
minf_fVar = [];
%
iterCount = 0;
fevalCount = 0;
startTime = time();
%
proglog_lastTime = time();
proglog_stopTime = time();
%
running_fevalCount = 0;
running_fTot = 0.0;
running_xtgTot = 0.0;
running_vecGTot = zeros(sizeX,1);
running_vecXTot = zeros(sizeX,1);
running_fSqTot = 0.0
running_xtgSqTot = 0.0;
running_vecGSqTot = zeros(sizeX,1);
running_vecXSqTot = zeros(sizeX,1);
%
record_matX = [];
record_matG = [];
record_rvecF = [];
record_rvecW = [];
record_indexToDrop = [];
record_matV = [];
if (prm.useQNJ)
	qnj_sMax = prm.qnj_sMaxInit;
	qnj_dMax = prm.qnj_dMaxInit;
endif
%
while ( 1 )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CONVENTIONAL GRADIENT DESCENT + STUFF
	%
	[ f, vecG ] = funchFG( vecX );
	fevalCount++;
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
	if ( running_fevalCount < prm.numFevalPerSuperPt )
		continue
	endif
	vecXHarvest = vecX;
	vecPHarvest = vecP;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SUPER-POINT ANALYSIS
	%
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
	elseif ( (superPt_f <= minf_f + 2.0*superPt_fVar + 2.0*minf_fVar) && (norm(superPt_vecG) < norm(best_vecG)) )
		newIsBest = true;
	endif
	if ( newIsBest )
		best_vecX = superPt_vecX;
		best_vecG = superPt_vecG;
		best_f = superPt_f;
		best_fVar = superPt_fVar;
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
	% PREPARE FOR NEXT ITER
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
	if ( ~prm.useQNJ )
		%
		vecXSeed = vecXHarvest;
		vecPSeed = vecPHarvest;
		vecX = vecXSeed;
		vecP = vecPSeed;
		continue;
	endif
	%
	error( "QNJ not yet supported." );
	%
	qnj_matV = [];
	record_indexToDrop = [];
	vecXSeed = vecXHarvest;
	vecPSeed = vecPHarvest;
	vecX = vecXSeed;
	vecP = vecPSeed;
endwhile
vecXF = best_vecX;
z__proglog0113;
z__fin0113;
