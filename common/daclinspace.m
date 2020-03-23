%  Function...
%    [ xVals, retCode, datOut ] = daclinspace( x0, x1, numValsRequested, funchY, prm=[] )
%  Overview...
%    Attempts to return a row vector with numValsRequested elements which produce
%     values of funcF which have an evenly spaced distance along the curve (DAC).
%    Compare to the built-in functions linspace() and logspace().
%  Input values...
%    x0, x1: The bounding values in the argument for funchY.
%    numValsRequested: The desired number of argument values, counting x0 and x1.
%    funchY: A function handle from R^1 to R^sizeY.
%    prm: Structure of parameters for the calculation.
%  Output values...
%    xVals: The output row vector of arguments to funchY.
%    retCode: A common return code, RETCODE__SUCCESS (0) on success.
%    datOut: Additional output data.
%  See source code for more information on prm and datOut.
%
%  This function is intended for use with postprocessing visualiation,
%   and may be inefficient and imprecise.
function [ rvecX, retCode, datOut ] = daclinspace( ...
  x0, x1, numValsRequested, funchY, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	commoninit;
	thisFile = "daclinspace";
	%
	assert( isrealscalar(x0) );
	assert( isrealscalar(x1) );
	assert( x0 ~= x1 );
	assert( isposintscalar(numValsRequested) );
	if ( 0 >= numValsRequested )
		rvecX = [];
		datOut.matY = [];
		return;
	elseif ( 1 == numValsRequested )
		rvecX = (x0+x1)/2.0;
		datOut.matY = funchY(rvecX);
		return;
	elseif ( 2 == numValsRequested )
		rvecX = [ x0, x1 ];
		datOut.matY(:,1) = funchY(x0);
		datOut.matY(:,2) = funchY(x1);
		return;
	end
	assert( 3 <= numValsRequested );
	%
	coeffMin = mygetfield(prm,"coeffMin",0.8);
	coeffMax = mygetfield(prm,"coeffMax",1.2);
	assert( isrealscalar(coeffMin) );
	assert( isrealscalar(coeffMax) );
	assert( 0.0 < coeffMin );
	assert( coeffMin < 1.0 );
	assert( 1.0 < coeffMax );
	%
	numIterLimit = mygetfield(prm,"numIterLimit",10);
	assert( isposintscalar(numIterLimit+1) );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	rvecX = linspace( x0, x1, numValsRequested );
	numIter = 0;
	while (1)
		numVals = size(rvecX,2);
		assert( isrealarray(rvecX,[1,numVals]) );
		assert( x0 <= rvecX );
		assert( x1 >= rvecX );
		rvecDeltaX = rvecX(2:end)-rvecX(1:end-1);
		assert( (x1-x0)*rvecDeltaX(:) > 0.0 );
		%
		if (mygetfield(prm,"funchYSupportsMultiArg",false))
			matY = funchY(rvecX);
		else
			clear matY;
			for n=1:size(rvecX,2)
				matY(:,n) = funchY(rvecX(n));
			end
		end
		sizeY = size(matY,1);
		assert( isrealarray(matY,[sizeY,numVals]) );
		%
		% Note that these deltaDAC values are simply linear estimates;
		% something like a spline would be more accurate.
		rvecDeltaDAC = sqrt(sum( (matY(:,2:end)-matY(:,1:end-1)).^2, 1 ));
		assert( rvecDeltaDAC(:) > 0.0 );
		rvecDAC = [ 0.0, cumsum(rvecDeltaDAC) ];
		fullDAC = rvecDAC(end);
		assert( fullDAC > 0.0 );
		%
		rvecDACDesired = linspace( 0.0, fullDAC, numValsRequested );
		deltaDACDesired = fullDAC / (numValsRequested-1.0);
		%
		if ( max(rvecDeltaDAC) <= deltaDACDesired * coeffMax ...
		  && min(rvecDeltaDAC) >= deltaDACDesired * coeffMin )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged in %d iterations ( %g ~ %g ).", ...
			  numIter, ...
			  min(rvecDeltaDAC), ...
			  max(rvecDeltaDAC) ) );
		  	retCode = RETCODE__SUCCESS;
			datOut.matY = matY;
		  	return;
		end
		if ( numIter >= numIterLimit )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Failed to converge after %d iterations ( %g ~ %g ).", ...
			  numIter, ...
			  min(rvecDeltaDAC), ...
			  max(rvecDeltaDAC) ) );
			retCode = RETCODE__IMPOSED_STOP;
			datOut.matY = matY;
			return;
		end
		%
		%
		rvecX_old = rvecX;
		% We're going to calculate a new rvecX.
		% This means we're going to throw out all of the matY,
		% and re-evaluate everything; it's inefficient, but easy to code.
		% We'll use simple linear interpolation between rvecDAC and
		% rvecDACDesired to get new values of rvecX.
		clear rvecN;
		n = 1;
		epsDAC = (eps^0.75)*fullDAC;
		for m=1:numVals
			while ( rvecDAC(n+1) + epsDAC < rvecDACDesired(m) );
				n++;
			end
			rvecN(m) = n;
		end
		rvecS = rvecDACDesired - rvecDAC(rvecN);
		rvecA = rvecS .* rvecDeltaX(rvecN) ./ rvecDeltaDAC(rvecN);
		rvecX = rvecX_old(rvecN) + rvecA;
		rvecX(1) = x0;
		rvecX(end) = x1;
		%
		%
		numIter++;
	end
end

%!test
%!	test_daclinspace
