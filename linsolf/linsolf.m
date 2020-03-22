%  Function...
%    function [ vecX, retCode, datOut ] ...
%      = linsolf( funchMatAProd, vecB, prm=[], datIn=[] )
%  Overview...
%    Finds an (inexact) solution to the equation matA * vecX = vecB
%    (essentially) using GMRes. The matrix A must be square, but
%    does not need to be provided explicitly, making this suitable
%    for JFNK.
%  Input...
%    funchMatAProd: A function handle like @(vecDummy)( matA * vecDummy ).
%    vecB: A column vector of the appropriate size.
%    prm: Input parameters.
%    datIn: Input data.
%  Output...
%    vecX: The obtained solution vector.
%    retCode: A standard return code, 0 indicates success.
%    datOut: Output data.
%  See source code for more information on prm, datIn, and datOut.
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
	reportInterval = mygetfield( prm, "reportInterval", 3.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	% Problem description.
	sizeF = size(vecB,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecB,[sizeF,1]) );
	sizeX = sizeF;
	normB = sqrt(sum(vecB.^2));
	%
	% Stopping criteria.
	fracResTol = mygetfield( prm, "fracResTol", 0.1 ); % Success.
	stopsigCheckInterval = mygetfield( prm, "stopsigCheckInterval", 1.0 ); % Imposed stop.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	numIterLimit = mygetfield( prm, "numIterLimit", -1.0 ); % Imposed stop.
	assert( isrealscalar(fracResTol) );
	assert( isrealscalar(stopsigCheckInterval) );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(numIterLimit) );
	stopsigCheckTimePrev = startTime;
	%
	% Internal parameters.
	gsThresh0 = mygetfield( prm, "gsThresh0", eps^0.50 );
	gsThresh1 = mygetfield( prm, "gsThresh1", eps^0.25 );
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
	matV = zeros(sizeX,0);
	matW = zeros(sizeF,0);
	matH = zeros(0,0);
	matR = zeros(0,0);
	vecG = zeros(0,1);
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
	vecBeta = vecB / normB;
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
		if ( 0.0 <= stopsigCheckInterval )
		if ( time() > stopsigCheckTimePrev + stopsigCheckInterval )
			if ( stopsignalpresent() )
				msg_notify( verbLev, thisFile, __LINE__, ...
				  sprintf("Found stop signal file in working directory.") );
				retCode = RETCODE__IMPOSED_STOP;
				linsolf__finish;
				return;
			end
			stopsigCheckTimePrev = time();
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
			   fracRes )  );
			reportTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		% Get new basis vector.
		if ( 0 == numIter )
			vecV = vecBeta;
		else
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
		assert( isrealarray(vecW,[sizeF,1]) );
		matW(:,numIter+1) = vecW;
		vecH = matW' * vecW;
		matH(1:numIter+1,numIter+1) = vecH;
		matH(numIter+1,1:numIter+1) = vecH';
		normW = sqrt(vecH(numIter+1,1));
		if ( 0.0 >= normW )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Projected vector %d is zero.", numIter+1) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			linsolf__finish;
			return;
		end
		[ matR, cholInfo ] = cholinsert( matR, numIter+1, vecH );
		if ( 0 != cholInfo )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Cholesky factorization %d failed.",numIter+1) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			linsolf__finish;
			return;
		end
		%
		vecG(numIter+1,1) = vecW' * vecBeta;
		vecZ = matR' \ vecG;
		fracRes = sqrt(max([ 0.0, 1.0-sum(vecZ.^2) ]));
		%
		% Prepare next iteration.
		numIter++;
		assert( sizeX == sizeF );
		assert( 0.0 < normW );
		vecV = matW(:,numIter) / normW;
	end
	%
	% Speed note:
	%  As long as we assume matV is always the Krylov subspace,
	%  we could get a speed-up by (essentially) reusing factorization
	%  information for V and W.
return;
end

%!test
%!	test_linsolf
