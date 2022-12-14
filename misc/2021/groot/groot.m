function [ vecX, retCode, datOut ] = groot( funchF, vecX0, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "groot";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	% Problem description.
	sizeX = size(vecX0,1);
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	funchOmegaOfF = @(f)( 0.5*(sum(f.^2,1)) );
	omega0 = funchOmegaOfF(vecF0);
	%
	% Stopping criteria.
	omegaTol = mygetfield( prm, "omegaTol", 1E-24 ); % Success.
	stopsigCheckInterval = mygetfield( prm, "stopsigCheckInterval", 1.0 ); % Imposed stop.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	numIterLimit = mygetfield( prm, "numIterLimit", -1.0 ); % Imposed stop.
	assert( isrealscalar(omegaTol) );
	assert( isrealscalar(stopsigCheckInterval) );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(numIterLimit) );
	stopsigCheckTimePrev = startTime;
	%
	% Internal parameters.
	%stepType_default = STEPTYPE__NEWTON;
	stepType_default = STEPTYPE__LEVCURVE;
	%stepType_default = STEPTYPE__GRADDIR;
	%stepType_default = STEPTYPE__GRADCURVE;
	stepTypeList = mygetfield( prm, "stepTypeList", [stepType_default] );
	funchJ = prm.funchJ;
	btIterLimit = 10;
	fallThresh = 1.0e-4;
	numFigs = 0;
	searchIterLimit = 10;
	numSearchVals = 21;
	searchOmegaTol = 1.001;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PREPARATIONAL WORK
	%
	% Initialize iterates.
	numIter = 0;
	vecX = vecX0;
	vecF = vecF0;
	omega = omega0;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	while (true)
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CHECK STOP
		%
		% Check "success" condition first.
		if ( 0.0 <= omegaTol )
		if ( omega <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", omega, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			groot__finish;
			return;
		end
		end
		%
		% Check "imposed stop" conditions.
		if ( 0 <= numIterLimit )
		if ( numIter >= numIterLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached numIterLimit (%d).",numIterLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot__finish;
			return;
		end
		end
		%
		if ( 0.0 <= exeTimeLimit )
		if ( time() >= startTime + exeTimeLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached exeTimeLimit (%g).",exeTimeLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot__finish;
			return;
		end
		end
		%
		if ( 0.0 <= stopsigCheckInterval )
		if ( time() > stopsigCheckTimePrev + stopsigCheckInterval )
			if ( stopsignalpresent() )
				msg_notify( verbLev, thisFile, __LINE__, ...
				  sprintf("Found stop signal file in working directory.") );
				retCode = RETCODE__IMPOSED_STOP;
				groot__finish;
				return;
			end
			stopsigCheckTimePrev = time();
		end
		end
		%
		if ( 1 <= numIter )
		if ( omega >= omegaPrev )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Failed to decrease omega (%g >= %g).",omega,omegaPrev) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			groot__finish;
			return;
		end
		if ( abs(omegaPrev-omega) <= abs(omegaPrev*fallThresh) )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Failed to decrease omega sufficiently (%g, %g).",omega,omegaPrev) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			groot__finish;
			return;
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% REPORT
		%
		if ( 0.0 <= reportInterval )
		if ( time() > reportTimePrev + reportInterval )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "  %4d  %6.2f  %8.2e", ...
			   numIter, ...
			   time()-startTime, ...
			   omega )  );
			reportTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		matJ = funchJ( vecX );
		vecG = matJ' * vecF;
		matH = matJ' * matJ;
		matI = eye(sizeF,sizeX);
		i0 = min([ numIter+1, max(size(stepTypeList)) ]);
		stepType = stepTypeList(i0);
		%
		datOut.iterDat(numIter+1).numIter = numIter;
		datOut.iterDat(numIter+1).vecX = vecX;
		datOut.iterDat(numIter+1).vecF = vecF;
		datOut.iterDat(numIter+1).omega = omega;
		datOut.iterDat(numIter+1).matJ = matJ;
		datOut.iterDat(numIter+1).vecG = vecG;
		datOut.iterDat(numIter+1).matH = matH;
		datOut.iterDat(numIter+1).matI = matI;
		datOut.iterDat(numIter+1).stepType = stepType;
		%
		%
		switch( stepType )
		case {STEPTYPE__NEWTON}
			vecT = -matH \ vecG;
			funchDelta = @(nu)( nu*vecT);
		case {STEPTYPE__GRADDIR}
			vecT = -vecG;
			t0 = vecT' * matH * vecT;
			if ( t0 > 0.0 )
				s0 = -(vecT'*vecG)/t0;
			else
				s0 = 0.0;
			end
			vecT *= s0;
			funchDelta = @(nu)( nu*vecT );
		case {STEPTYPE__LEVCURVE}
			hScl = max(diag(matH));
			assert( 0.0 < hScl );
			matA = (matH/hScl) - matI;
			vecT = -(vecG/hScl);
			funchDelta = @(nu)( nu*( (matI+(nu*matA))\vecT ) );
		case {STEPTYPE__GRADCURVE}
			[ matPsi, matLambda ] = eig( matH );
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeX^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeX^3)*(eps^0.75) );
			vecPsiTN = matPsi'*(-matH\vecG);
			lambdaMin = min(diag(matLambda));
			matSigma = matLambda / lambdaMin;
			funchDelta = @(nu)( ...
			  matPsi * ( vecPsiTN - (diag(nu.^diag(matSigma))*vecPsiTN) ) );
		otherwise
			error(sprintf("Unsupported value of stepTypeList(%d) (%d).", ...
			  i0, stepTypeList(i0) ));
		end
		%
		%
		funchDeltaNorm = @(nu)( sqrt(sum((funchDelta(nu)).^2)) );
		searchIter = 0;
		nuMin = 0.0;
		nuMax = 1.0;
		while (true)
			%echo__nuRange = [ nuMin, nuMax ]
			rvecNu = linspace( nuMin, nuMax, numSearchVals );
			clear matDelta;
			for n=1:numSearchVals
				matDelta(:,n) = funchDelta(rvecNu(n));
			end
			rvecDeltaNorm = sqrt(sum(matDelta.^2,1));
			matF = funchF(repmat(vecX,[1,numSearchVals])+matDelta);
			rvecOmega = funchOmegaOfF(matF);
			%for n=1:numSearchVals
			%	rvecDeltaNorm(n) = sqrt(sum(matDelta(:,n).^2,1));
			%	vecF = funchF(vecX+matDelta(:,n));
			%	rvecOmega(n) = funchOmegaOfF(vecF);
			%end
			%
			if (0)
			numFigs++; figure(numFigs);
			plot( rvecNu, rvecDeltaNorm, 'o-' );
			title(sprintf( "deltaNorm v nu %d x %d", numIter, searchIter ));
			grid on;
			%
			numFigs++; figure(numFigs);
			plot( rvecDeltaNorm, rvecOmega, 'o-' );
			title(sprintf( "omega v deltaNorm %d x %d", numIter, searchIter ));
			grid on;
			end
			%
			[ omegaBest, nOfBest ] = min([ rvecOmega ]);
			%echo__omegaBest = omegaBest
			nuOfBest = rvecNu(nOfBest);
			searchIter++;
			if ( searchIterLimit <= searchIter )
				break;
			elseif ( max(rvecOmega) < searchOmegaTol*omegaBest )
				break;
			end
			if (1==nOfBest)
				nuMax = rvecNu(3);
			elseif (numSearchVals==nOfBest)
				nuMin = rvecNu(numSearchVals-2);
			else
				nuMin = rvecNu(nOfBest-1);
				nuMax = rvecNu(nOfBest+1);
			end
		end
		%
		%echo__nuOfBest = nuOfBest
		vecDelta = funchDelta(nuOfBest);
		vecXTrial = vecX + vecDelta;
		vecFTrial = funchF( vecXTrial );
		omegaTrial = funchOmegaOfF( vecFTrial );
		%
		clear vecT;
		clear funchDelta;
		clear matA;
		clear hScl;
		vecXPrev = vecX;
		vecFPrev = vecF;
		omegaPrev = omega;
		%
		numIter++;
		vecX = vecXTrial;
		vecF = vecFTrial;
		omega = omegaTrial;
		datOut.iterDat(numIter).vecXTrial = vecXTrial;
		datOut.iterDat(numIter).vecFTrial = vecFTrial;
		datOut.iterDat(numIter).omegaTrial = omegaTrial;
	end
	%
return;
end

%!test
%!	test_groot
