% function [ xOfCand, meritOfCand, datOut ] = extFinder( xVals, fVals, prm = [], datIn = [] );
	prm = [];
	datIn = [];
	%
	%plot( xVals, fVals, 'o-' );
	%grid on;
	%
	extFinder__init;
	thisFile = "extFinder";
	if (haveCand)
		thisFile = [ "RETURN from " thisFile ];
		return;
	end
	%
	%
	bigS0 = mygetfield( prm, "bigS0", [] );
	bigP0 = mygetfield( prm, "bigP0", [] );
	if ( isempty(bigS0) || isempty(bigP0) )
		prm_dlogfdf = mygetfield( prm, "prm_dlogfdf", [] );
		[ bigX_temp, bigP_temp ] = extFinder_dlogdf( xVals, fVals, prm_dlogfdf );
		if ( isempty(bigP0) )
			bigP0 = bigP_temp;
		end
		assert( isrealscalar(bigP0) );
		if ( isempty(bigS0) )
			bigP = bigP0;
			bigS = bigX_temp;
			[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, fVals, nOfPtWiseMin );
			xExt = bigS - bigDelta*mypow( bigB/(bigA*bigP), 1.0/(bigP-1.0) )
			assert( xVals(nOfPtWiseMin-1) < xExt ); % What to do if this fails?!?!
			assert( xExt < xVals(nOfPtWiseMin+1) ); % What to do if this fails?!?!
			bigS0 = bigX_temp;
		end
		assert( isrealscalar(bigS0) );
		clear bigX_temp;
		clear bigP_temp;
		clear prm_dlogfdf;
		%
		%
		echo__bigS0 = bigS0
		echo__bigP0 = bigP0
	end
	assert( isrealscalar(bigS0) );
	assert( isrealscalar(bigP0) );
	bigS = bigS0;
	bigP = bigP0;
	[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, fVals, nOfPtWiseMin );
	%
	%
	funchGModel = @(x)( bigG0 + bigG1*( ...
	  bigA*abs((x-bigS)/bigDelta).^bigP + bigB*(x-bigS)/bigDelta + bigC ) );
	plot( ...
	  xVals, gVals, 'o-', ...
	  xVals, funchGModel(xVals), 'x-' );
	grid on;
	%
	%
return;
