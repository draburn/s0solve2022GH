clear;
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
	return;
endif
stopsig_lastTime = time();
[ funchFG, vecX0, prm, initDat ] = z__init0112();
sizeX = size(vecX0,1);
%
vecX = vecX0;
vecP = zeros(sizeX,1);
vecXSeed = vecX;
vecPSeed = vecP;
%
best_vecX = vecX0;
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
while ( 1 )
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
	%
	msg( __FILE__, __LINE__, " Yo! This deserves more consideration..." );
	msg( __FILE__, __LINE__, sprintf( "  f = %0.3e +/- %0.3e (vs %0.3e)", superPt_fAvg, sqrt(superPt_fSqVar), superPt_f ) );
	msg( __FILE__, __LINE__, sprintf( "  x = %0.3e +/- %0.3e", norm(superPt_vecX), sqrt(sum(superPt_vecXSqVar)) ) );
	msg( __FILE__, __LINE__, sprintf( "  g = %0.3e +/- %0.3e", norm(superPt_vecG), sqrt(sum(superPt_vecGSqVar)) ) );
	msg( __FILE__, __LINE__, sprintf( "  xtg = %0.3e +/- %0.3e", norm(superPt_xtgAvg), sqrt(sum(superPt_xtgSqVar)) ) );
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
	best_vecX = superPt_vecX; % This is a big assumption, but not entirely unreasonable.
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
	%
	if ( time() > proglog_lastTime + prm.progressReportInterval )
		z__proglog0112;
	 	proglog_lastTime = time();
	endif
	%
	vecXSeed = vecXHarvest;
	vecPSeed = vecPHarvest;
	vecX = vecXSeed;
	vecP = vecPSeed;
endwhile
vecXF = best_vecX;
z__proglog0112;
z__fin0112;
