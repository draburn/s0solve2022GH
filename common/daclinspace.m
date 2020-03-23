function [ rvecX, retCode, datOut ] = daclinspace( ...
  x0, x1, numValsRequested, funchF, prm=[], datIn )
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
		return;
	elseif ( 1 == numValsRequested )
		rvecX = (x0+x1)/2.0;
		return;
	elseif ( 2 == numValsRequested )
		rvecX = [ x0, x1 ];
		return;
	end
	assert( 3 <= numValsRequested );
	%
	tolMin = mygetfield(prm,"tolMin",0.6);
	tolMax = mygetfield(prm,"tolMax",1.4);
	assert( isrealscalar(tolMin) );
	assert( isrealscalar(tolMax) );
	assert( 0.0 < tolMin );
	assert( tolMin < 1.0 );
	assert( 1.0 < tolMax );
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
		if (mygetfield(prm,"funchFSupportsMultiArg",false))
			matF = funchF(rvecX);
		else
			clear matF;
			for n=1:size(rvecX,2)
				matF(:,n) = funchF(rvecX(n));
			end
		end
		sizeF = size(matF,1);
		assert( isrealarray(matF,[sizeF,numVals]) );
		%
		rvecDeltaDAC = sqrt(sum( (matF(:,2:end)-matF(:,1:end-1)).^2, 1 ));
		assert( rvecDeltaDAC(:) > 0.0 );
		rvecDAC = [ 0.0, cumsum(rvecDeltaDAC) ];
		fullDAC = rvecDAC(end);
		assert( fullDAC > 0.0 );
		%
		rvecDACDesired = linspace( 0.0, fullDAC, numValsRequested );
		%numFigs++; figure(numFigs);
		%plot( rvecDACDesired, rvecDAC, 'o-' );
		%grid on;
		%
		if ( min(rvecDeltaDAC) >= tolMin * fullDAC / numValsRequested ...
		  && max(rvecDeltaDAC) <= tolMax * fullDAC / numValsRequested )
		  	retCode = RETCODE__SUCCESS;
		  	echo__numIter = numIter
		  	return;
		end
		if ( numIter >= numIterLimit )
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		% There has got to be a better way to do this...
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
		%
		rvecX_old = rvecX;
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
