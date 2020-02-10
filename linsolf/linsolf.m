function [ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "linsolf";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 1.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - reportInterval - 1.0;
	%
	% Problem description.
	sizeF = size(vecB,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecB,[sizeF,1]) );
	sizeX = sizeF;
	normB = sqrt(sum(vecB.^2));
	%
	% Stopping criteria.
	fracResTol = mygetfield( prm, "fracResTol", 1.0e-1 ); % Success.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	numIterLimit = mygetfield( prm, "numIterLimit", -1.0 ); % Imposed stop.
	assert( isrealscalar(fracResTol) );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(numIterLimit) );
	%
	% Internal parameters.
	gsThresh0 = mygetfield( prm, "gsThresh0", eps^0.75 );
	gsThresh1 = mygetfield( prm, "gsThresh1", 0.5 );
	assert( isrealscalar(gsThresh0) );
	assert( 0.0 < gsThresh0 );
	assert( 1.0 > gsThresh0 );
	assert( isrealscalar(gsThresh1) );
	assert( 0.0 < gsThresh1 );
	assert( 1.0 > gsThresh1 );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PREPARATIONAL WORK
	%
	% Initialize iterates.
	numIter = 0;
	matU = zeros(sizeX,0);
	matV = zeros(sizeX,0);
	matW = zeros(sizeF,0);
	vecX = zeros(sizeX,1);
	%
	% Handle the "start is solution" case.
	if ( 0.0 >= normB )
		fracRes = 0.0;
		msg_notify( verbLev, thisFile, __LINE__, "Initial residual is already zero." );
		retCode = RETCODE__SUCCESS;
		linsolf__finish;
		return;
	end
	fracRes = 1.0;
	assert( sizeX == sizeF );
	matU(:,1) = vecB / normB;
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
		if ( 0.0 <= fracResTol )
		if ( fracRes <= fracResTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", fracRes, fracResTol ) );
			retCode = RETCODE__SUCCESS;
			linsolf__finish;
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
			linsolf__finish;
			return;
		end
		end
		%
		if ( 0.0 <= exeTimeLimit )
		if ( time() >= startTime + exeTimeLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached exeTimeLimit (%g).",exeTimeLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			linsolf__finish;
			return;
		end
		end
		%
		if ( stopsignalpresent() )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Found stop signal file in working directory.") );
			retCode = RETCODE__IMPOSED_STOP;
			linsolf__finish;
			return;
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% REPORT
		%
		if ( verbLev >= VERBLEV__PROGRESS )
		if ( time() > reportTimePrev + reportInterval )
			msg( thisFile, __LINE__, sprintf( ...
			  "  %4d  %6.2f  %8.2e", ...
			   numIter, ...
			   time()-startTime, ...
			   fracRes )  );
			reportTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		% Get new basis vector.
		vecV = matU(:,numIter+1);
		assert( isrealarray(vecV,[sizeX,1]) );
		assert( abs(1.0-sum(vecV.^2)) < eps^0.75 );
		if ( 1 <= numIter )
			vecV -= matV * ( matV' * vecV );
			normV = sqrt(sum(vecV.^2));
			if ( normV <= gsThresh0 )
				msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
				  "Basis vector %d was destroyed during first pass (%e < %e).", ...
				  numIter+1, normV, gsThresh0) );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				linsolf__finish;
				return;
			end
			vecV /= normV;
			%
			vecV -= matV * ( matV' * vecV );
			normV = sqrt(sum(vecV.^2));
			if ( normV <= gsThresh1 )
				msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
				  "Basis vector %d is very small after second pass (%e < %e).", ...
				  numIter+1, normV, gsThresh1) );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				linsolf__finish;
				return;
			end
			vecV /= normV;
		end
		matV(:,numIter+1) = vecV;
		%
		% Get next projected vector...
		vecW = funchMatAProd( vecV );
		matW(:,numIter+1) = vecW;
		normW = sqrt(sum(vecW.^2));
		if ( 0.0 >= normW )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Projected vector %d is zero.", numIter+1) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			linsolf__finish;
			return;
		end
		%
		% Get new solution.
		% This can be done way, way more efficiently.
		%  In particular (see s0solve20200104)...
		%   1. Use a Cholesky factorization of W^T W;
		%   2. Use of cholinsert() to update the factorization; and,
		%   3. Use a more direct calculation of resFrac, like...
		%    vecG(numIter+1,1) = vecW' * vecBeta;
		%    vecZ = matR' \ vecG;
		%    resFrac = sqrt(max([ 0.0, 1.0 - sum(vecZ.^2) ]));
		%  Note that, when V is *always* the Krylov subspace,
		%   there's an additional opportunity for speed-up from
		%   (essentially) not having to orthonormalize both V and W.
		%   But, I prefer to allow non-Krylov flexibility.
		% But, premature optimization is the root of all evil.
		vecY = matW \ vecB;
		vecX = matV * vecY;
		vecRho = (matW * vecY) - vecB;
		res = sqrt(sum(vecRho.^2));
		fracRes = res / normB;
		%
		% Prepare next iteration.
		numIter++;
		assert( sizeX == sizeF );
		assert( 0.0 < normW );
		matU(:,numIter+1) = matW(:,numIter) / normW;
	end
	%
return;
end
