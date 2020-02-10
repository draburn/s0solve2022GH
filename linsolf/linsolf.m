function [ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "linsolf";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	datOut.funchMatAProd = funchMatAProd;
	datOut.vecB = vecB;
	datOut.prm = prm;
	datOut.thisFile = thisFile;
	datOut.startTime = startTime;
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
	datOut.verbLev = verbLev;
	datOut.reportInterval = reportInterval;
	reportTimePrev = startTime - abs(reportInterval) - 1.0;
	%
	% Problem description.
	sizeF = size(vecB,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecB,[sizeF,1]) );
	sizeX = sizeF;
	normB = sqrt(sum(vecB.^2));
	datOut.sizeF = sizeF;
	datOut.sizeX = sizeX;
	datOut.normB = normB;
	%
	% Stopping criteria.
	fracResTol = mygetfield( prm, "fracResTol", 1.0e-2 ); % Success.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	numIterLimit = mygetfield( prm, "numIterLimit", 100 ); % Imposed stop.
	assert( isrealscalar(fracResTol) );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(numIterLimit) );
	datOut.fracResTol = fracResTol;
	datOut.exeTimeLimit = exeTimeLimit;
	datOut.numIterLimit = numIterLimit;
	%
	% Internal parameters.
	orthonormThresh0 = mygetfield( prm, "orthonormThresh0", eps^0.75 );
	orthonormThresh1 = mygetfield( prm, "orthonormThresh1", 0.5 );
	assert( isrealscalar(orthonormThresh0) );
	assert( 0.0 < orthonormThresh0 );
	assert( 1.0 > orthonormThresh0 );
	assert( isrealscalar(orthonormThresh1) );
	assert( 0.0 < orthonormThresh1 );
	assert( 1.0 > orthonormThresh1 );
	datOut.orthonormThresh0 = orthonormThresh0;
	datOut.orthonormThresh1 = orthonormThresh1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PREPARATIONAL WORK
	%
	% Initialize iterates.
	numIter = 0;
	sizeK = 0;
	matU = zeros(sizeX,0);
	matV = zeros(sizeX,0);
	matW = zeros(sizeF,0);
	vecX = zeros(sizeX,1);
	%
	% Update datOut.
	datOut.numIter = numIter;
	datOut.sizeK = sizeK;
	datOut.matU = matU;
	datOut.matV = matV;
	datOut.matW = matW;
	datOut.vecX = vecX;
	%
	% Handle the corner case of vecB = 0.
	if ( 0.0 >= normB )
		msg_notify( verbLev, thisFile, __LINE__, "Initial residual is already zero." );
		retCode = RETCODE__SUCCESS;
		return;
	end
	fracRes = 1.0;
	vecBeta = vecB / normB;
	assert( sizeX == sizeF );
	matU(:,1) = vecBeta;
	datOut.fracRes = fracRes;
	datOut.vecBeta = vecBeta;
	datOut.matU = matU;
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
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Final:  %4d  %6.2f  %8.2e  %s", ...
			   numIter, ...
			   time()-startTime, ...
			   fracRes, ...
			   retcode2str(retCode) )  );
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
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Final:  %4d  %6.2f  %8.2e  %s", ...
			   numIter, ...
			   time()-startTime, ...
			   fracRes, ...
			   retcode2str(retCode) )  );
			return;
		end
		end
		%
		if ( 0.0 <= exeTimeLimit )
		if ( time() >= startTime + exeTimeLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached exeTimeLimit (%g).",exeTimeLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Final:  %4d  %6.2f  %8.2e  %s", ...
			   numIter, ...
			   time()-startTime, ...
			   fracRes, ...
			   retcode2str(retCode) )  );
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
		msg_error( verbLev, thisFile, __LINE__, "ERROR: NOT IMPLEMENTED!" );
		return;
		%
		numIter++;
		sizeK++;
	end
	%
return;
end
