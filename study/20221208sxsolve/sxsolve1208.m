function [ vecXBest, retCode, datOut ] = sxsolve1208( funchFG, vecXInit, prm=[] )
	mydefs;
	startTime = time();
	vecXBest = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	sizeX = size(vecXInit,1);
	assert( isposintscalar(sizeX) );
	assert( isrealarray(vecXInit,[sizeX,1]) );
	prm = __init( funchFG, vecXInit, prm );
	if ( prm.stopSignalCheckInterval >= 0.0 )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		stopSignalCheckedTime = time();
	endif
	%
	[ f0, vecG0 ] = funchFG( vecXInit );
	fevalCount = 1;
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	matX = vecXInit;
	rvecF = f0;
	matG = vecG0;
	iterCount = 0;
	%indexOfBest = 1;
	%
	vecXBestPrev = vecXInit;
	fBestPrev = f0;
	vecGBestPrev = vecG0;
	vecXBest = vecXInit;
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
			matD = matX - vecXBest;
			matV = utorthdrop( matD, prm.dropThresh );
			sizeK = size(matV,2);
			matVTD = matV' * matD;
			matVTG = matV' * matG;
			[ fFit, vecVTGFit, matVTHVFit ] = hessfit( matVTD, rvecF, matVTG );
			if ( size(matV,2) >= sizeX )
				echo__matHFit = matV * matVTHVFit * (matV')
			endif
			%matD
			%matV
			%matVTD
			%matVTG
			%matG
			matG_before = matG;
			check_matG_before = matV*( vecVTGFit + matVTHVFit*matVTD );
			%rd_matG = reldiff( check_matG_before, matG )
			%
			vecGPerp = vecGBest - matV * ( matV'*vecGBest );
			normGPerp = norm(vecGPerp);
			if ( normGPerp > prm.dropThresh*norm(vecGBest) )
				vecGHat = vecGPerp/normGPerp;
				matV = [ matV, vecGHat ];
				matVTD = matV' * matD;
				matVTG = matV' * matG;
				%vecGPerp
				%matV
				%matD
				%matG
				%matVTD
				%matVTG
				%vecVTGFit
				matVTHVFit
				eigH_before = eig(matVTHVFit)
				[ vecVTGFit, matVTHVFit ] = __appendFit( matVTD, matVTG, vecVTGFit, matVTHVFit );
				%vecVTGFit
				matVTHVFit
				eigH_after = eig(matVTHVFit)
				%
				check_matG_before
				matG_before
				matG
				%check_matG = matV*( vecVTGFit + matVTHVFit*matVTD )
				%rd_matG = reldiff( check_matG, matG )
				error( "Begone!" );
			endif
			%
			%matB = diag(max( abs(matVTD), [], 2 ));
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
			stepFactor *= prm.btCoeff;
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


function prm = __init( funchFG, vecXInit, prmIn )
	mydefs;
	sizeX = size(vecXInit,1);
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
	prm.gradStepCoeff = 0.1;
	prm.dropThresh = eps^0.7;
	prm.btCoeff = 0.1;
	%
	prm.numPtPerSuperPt = 100;
	prm.momentumCoeff = 0.9;
	prm.learningRate = 0.01;
	%
	prm = overwritefields( prm, prmIn );
	assert( prm.smallFallTrgt > 0.0 );
return;
endfunction


function [ fevalCount, vecXSuperPt, fSuperPt, vecGSuperPt, weightSuperPt, vecXFin, vecPFin ] = __evalSuperPt( funchFG, vecXInit, vecPInit, prm )
	fevalCount = 0;
	vecXSuperPt = vecXInit;
	xtgSuperPt = 0.0;
	fSuperPt = 0.0;
	vecGSuperPt = zeros(size(vecXInit));
	weightSuperPt = 0.0;
	vecX = vecXInit;
	vecP = vecPInit;
	for p = 1 : prm.numPtPerSuperPt
		[ f, vecG ] = funchFG( vecX );
		fevalCount++;
		vecXSuprtPt += vecX;
		fSuperPt += f;
		xtgSuperPt += vecX'*vecG;
		vecGSuperPt += vecG;
		weightSuperPt += 1.0;
		vecP = prm.momentumCoeff * vecP + prm.learningRate * vecG;
		vecX += vecP;
	endfor
	vecXSuperPt /= weightSuperPt;
	vecGSuperPt /= weightSuperPt;
	xtgSuperPt /= weightSuperPt;
	fSuperPt = (fSuperPt/prm.numPtPerSuperPt) - 0.5 *( xtgSuperPt - (vecXSuperPt'*vecGSuperPt) );
	vecXFin = vecX;
	vecPFin = vecP;
return;
endfunction


function [ vecVTGFit, matVTHVFit ] = __appendFit( matVTD, matVTG, vecVTGFit, matVTHVFit )
	sizeK_before = size(matVTHVFit,1)
	%
	appendTo_vecVTGFit = matVTG(end);
	%
	[ foo, indexAnchor ] = min( sum(matVTD.^2, 1) );
	assert( foo < eps*max( sum(matVTD.^2, 1) ) );
	%
	%matY = [ matVTD(:,1:indexAnchor-1), matVTD(:,indexAnchor+1:end) ];
	%matY = [ matY(1:indexAnchor-1,:); matY(indexAnchor+1:end,:) ];
	%matO = [ matVTG(:,1:indexAnchor-1), matVTG(:,indexAnchor+1:end) ];
	%vecO = matO(indexAnchor,:)';
	matY = [ matVTD(:,1:indexAnchor-1), matVTD(:,indexAnchor+1:end) ];
	matY = matY(1:end-1,:);
	matO = [ matVTG(:,1:indexAnchor-1), matVTG(:,indexAnchor+1:end) ];
	vecO = matO(indexAnchor,:)';
	%
	vec_appendTo_matVTHVFit = matY \ (vecO - appendTo_vecVTGFit);
	%vec_appendTo_matVTHVFit = zeros( sizeK_before, 1 );
	%
	h = sqrt( sumsq(vec_appendTo_matVTHVFit) + sumsq(diag(matVTHVFit))/sizeK_before );
	%h =  sqrt(sumsq(vec_appendTo_matVTHVFit)) + sqrt(sumsq(diag(matVTHVFit)));
	vecVTGFit = [ vecVTGFit; appendTo_vecVTGFit ];
	matVTHVFit = [ matVTHVFit, vec_appendTo_matVTHVFit; vec_appendTo_matVTHVFit', h ];
return;
endfunction
