%  Function...
%    [ xVals, retCode, datOut ] = flinspace( x0, x1, numValsRequested, funchF, prm=[] )
%  Overview...
%    Attempts to return a row vector with numValsRequested elements which produce
%     evenly spaced elements when used as arguments to funchF.
%    Compare to the built-in functions linspace() and logspace().
%    If is assumed that funchF monotonic between x0 and x1;
%     if this is not the case, results may be strange.
%    It is preferred for x0 to be less than x1 and for funchF to be strictly
%     increasing, but these properties shouldn't have to be satisfied.
%  Input values...
%    x0, x1: The bounding values in the argument for funchF.
%    numValsRequested: The desired number of argument values, counting x0 and x1.
%    funchF: A function handle for a monotonic function.
%    prm: Structure of parameters for the calculation.
%  Output values...
%    xVals: The output row vector of arguments to funchF.
%    retCode: A common return code, RETCODE__SUCCESS (0) on success.
%    datOut: Additional output data.
%  See source code for more information on prm and datOut.
function [ xVals, retCode, datOut ] = flinspace( x0, x1, numValsRequested, funchF, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "flinspace";
	retCode = RETCODE__NOT_SET;
	datOut = [];
	startTime = time();
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	reportInterval = mygetfield( prm, "reportInterval", 3.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PROBLEM-SPECIFIC INIT.
	%
	assert( isrealscalar(x0) );
	assert( isrealscalar(x1) );
	assert( x0 ~= x1 );
	if ( x0 > x1 )
		temp0 = x0;
		x0 = x1;
		x1 = temp0;
		flippedX = true;
	else
		flippedX = false;
	end
	%
	assert( isrealscalar(numValsRequested) );
	assert( 0 <= numValsRequested );
	if ( 0 == numValsRequested )
		xVals = [];
		return;
	elseif ( 1 == numValsRequested )
		xVals = (x0+x1)/2.0;
		return;
	elseif ( 2 == numValsRequested )
		xVals = [ x0, x1 ];
		return;
	end
	assert( 3 <= numValsRequested );
	%
	f0 = funchF(x0);
	f1 = funchF(x1);
	assert( isrealscalar(f0) );
	assert( isrealscalar(f1) );
	assert( f0 ~= f1 );
	if ( f0 > f1 )
		funchF_input = funchF;
		funchF = @(dummy)( -funchF_input(dummy) );
		f0 = -f0;
		f1 = -f1;
		flippedF = true;
	else
		flippedF = false;
	end
	%
	funchDFDX = mygetfield( prm, "funchDFDX", [] );
	if ( isempty(funchDFDX) )
		epsDFDX = mygetfield( prm, "epsDFDX", abs(x1-x0)*(eps^0.75) );
		assert( isrealscalar(epsDFDX) );
		assert( epsDFDX > 0.0 );
		funchDFDX = @(dummy)( ...
		  ( funchF(min([x1,dummy+epsDFDX])) - funchF(max([x0,dummy-epsDFDX])) ) ...
		 /( min([x1,dummy+epsDFDX]) - max([x0,dummy-epsDFDX]) )  );
	else
		if (flippedF)
			funchDFDX_input = funchDFDX;
			funchDFDX = @(dummy)(-funchDFDX_input(dummy));
		end
	end
	%
	assert( x0 < x1 );
	assert( f0 < f1 );
	%
	deltaFTrgtOriginal = (f1-f0)/(numValsRequested-1.0);
	deltaFTrgtMinCoeff = mygetfield( prm, "deltaFTrgtMinCoeff", 0.6 );
	deltaFTrgtMaxCoeff = mygetfield( prm, "deltaFTrgtMinCoeff", 2.0 );
	deltaFMinCoeff = mygetfield( prm, "deltaFTrgtMinCoeff", 0.4 );
	deltaFMaxCoeff = mygetfield( prm, "deltaFTrgtMinCoeff", 3.0 );
	fTolCoeff = mygetfield( prm, "fTolCoeff", 0.01 );
	assert( isrealscalar(deltaFTrgtMinCoeff) );
	assert( isrealscalar(deltaFTrgtMaxCoeff) );
	assert( isrealscalar(deltaFMinCoeff) );
	assert( isrealscalar(deltaFMaxCoeff) );
	assert( isrealscalar(fTolCoeff) );
	assert( 0.0 <= fTolCoeff );
	assert( fTolCoeff < deltaFMinCoeff );
	assert( deltaFMinCoeff < deltaFTrgtMinCoeff );
	assert( deltaFTrgtMinCoeff < deltaFTrgtMaxCoeff );
	assert( deltaFTrgtMaxCoeff < deltaFMaxCoeff );
	%
	deltaFTrgtMin = deltaFTrgtMinCoeff * deltaFTrgtOriginal;
	deltaFTrgtMax = deltaFTrgtMaxCoeff * deltaFTrgtOriginal;
	deltaFMin = deltaFMinCoeff * deltaFTrgtOriginal;
	deltaFMax = deltaFMaxCoeff * deltaFTrgtOriginal;
	fTol = fTolCoeff * deltaFTrgtOriginal;
	%
	subIterLimit = mygetfield( prm, "subIterLimit", 50 );
	linearityThresh = mygetfield( prm, "linearityThresh", 0.01 );
	dfdxThreshCoeff = mygetfield( prm, "dfdxThreshCoeff", sqrt(eps) );
	assert( isrealscalar(subIterLimit) );
	assert( 1 <= subIterLimit );
	assert( isrealscalar(linearityThresh) );
	assert( 0.0 < linearityThresh );
	assert( isrealscalar(dfdxThreshCoeff) );
	assert( 0.0 <= dfdxThreshCoeff );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	xVals(1) = x0;
	fVals(1) = f0;
	numVals = 1;
	doMainLoop = false;
	while (true)
		xPrev = xVals(numVals);
		fPrev = fVals(numVals);
		if ( fPrev + deltaFTrgtMax >= f1 )
			break;
		end
		numNewValsTrgtMax = floor( (f1-fPrev)/deltaFTrgtMin );
		numNewValsTrgtMin = floor( (f1-fPrev)/deltaFTrgtMax );
		numNewValsTrgtIdeal = numValsRequested-numVals;
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "%2d; %9.2e, %9.2e; %2d, %2d, %2d", ...
		  numVals, xPrev, fPrev, ...
		  numNewValsTrgtMin, numNewValsTrgtIdeal, numNewValsTrgtMax ) );
		if ( numNewValsTrgtMax <= 1 )
			break;
		end
		assert( numNewValsTrgtMax >= numNewValsTrgtMin );
		numNewValsTrgt = median([ ...
		  numNewValsTrgtMin, numValsRequested-numVals, numNewValsTrgtMax ]);
		fTrgt = fPrev + ((f1-fPrev)/numNewValsTrgt);
		assert( fTrgt <= fPrev + deltaFTrgtMax );
		fMin = fPrev + deltaFMin;
		fMax = min([ f1, fPrev + deltaFMax ]);
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "%9.2e %9.2e, %9.2e, %9.2e, %9.2e; %9.2e, %9.2e", ...
		  fPrev, fMin, fTrgt, fMax, f1, xPrev, x1 ) );
		assert( fPrev < fMin );
		assert( fMin < fTrgt );
		assert( fTrgt < fMax );
		assert( fMax <= f1 );
		%
		if (verbLev >= VERBLEV__PROGRESS)
		if (2<=numVals)
		if (time()>reportTimePrev+reportInterval)
			msg( thisFile, __LINE__, sprintf( ...
			  "Estimated time remaining: %fs.", ...
			  (numNewValsTrgt-1.0) * (time()-startTime)/(numVals-1.0) ) );
			reportTimePrev = time();
		end
		end
		end
		xLo = xVals(numVals);
		fLo = fVals(numVals);
		xHi = x1;
		fHi = f1;
		%
		subIter = 0;
		while (true)
			% model: f = c0 + c1*(x-xLo) + c2*(x-xLo)^2.
			c0 = fLo;
			dx = xHi - xLo;
			if ( abs(fLo-fTrgt) < abs(fHi-fTrgt) )
				dfdxLo = funchDFDX(xLo);
				if ( dfdxLo < 0.0 )
				if ( abs(dfdxLo) < (abs(fLo)+abs(fHi))*dfdxThreshCoeff )
					dfdxLo = 0.0;
				end
				end
				if ( dfdxLo < 0.0 )
					errMsg = sprintf( ...
					  "ERROR: Derivative has wrong sign (%g, %g, %g, %g).", ...
					  f0, f1, xLo, dfdxLo );
					error( errMsg );
				end
				assert( dfdxLo >= 0.0 )
				c1 = dfdxLo;
				c2 = ( fHi - fLo - (dfdxLo*dx) )/( dx^2 );
			else
				dfdxHi = funchDFDX(xHi);
				if ( dfdxHi < 0.0 )
				if ( abs(dfdxHi) < (abs(fLo)+abs(fHi))*dfdxThreshCoeff )
					dfdxHi = 0.0;
				end
				end
				if ( dfdxHi < 0.0 )
					errMsg = sprintf( ...
					  "ERROR: Derivative has wrong sign (%g, %g, %g, %g).", ...
					  f0, f1, xHi, dfdxHi );
					error( errMsg );
				end
				assert( dfdxHi >= 0.0 );
				c2 = ( dfdxHi - ((fHi-fLo)/dx) ) / dx;
				c1 = dfdxHi - (2.0*c2*dx);
			end
			msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
			  "  %2d, %9.2e, %9.2e; %9.2e, %9.2e; %9.2e, %9.2e, %9.2e", ...
			  subIter, fLo, fHi, xLo, xHi, c0, c1, c2 ) );
			if ( abs(c2*dx) < linearityThresh*abs(c1) )
				deltaX = ((fTrgt-fLo)*(xHi-xLo)/(fHi-fLo));
			elseif (c2>0.0)
				discr = (c1^2) + (4.0*c2*(fTrgt-c0));
				assert( 0.0 < discr );
				deltaX = ( -c1 + sqrt(discr) ) / (2.0*c2);
			else % c2<0.0
				discr = (c1^2) + (4.0*c2*(fTrgt-c0));
				assert( 0.0 < discr );
				deltaX = ( -c1 + sqrt(discr) ) / (2.0*c2);
			end
			msg_copious( verbLev, thisFile, __LINE__, sprintf("  %9.2e", deltaX) );
			assert( 0.0 < deltaX );
			assert( xLo + deltaX < xHi );
			xNew = xLo + deltaX;
			fNew = funchF(xNew);
			if ( (fNew<fLo) || (fNew>fHi) )
				errMsg = sprintf( ...
				  "ERROR: Function is not monotonic: ( %g, %g, %g, %g, %g, %g ).", ...
				  xLo, xNew, xHi, fLo, fNew, fHi );
				error( errMsg );
			end
			if ( fNew < fTrgt )
				fLo = fNew;
				xLo = xNew;
			else
				fHi = fNew;
				xHi = xNew;
			end
			%
			% Do a bisection.
			xNew = (xLo+xHi)/2.0;
			fNew = funchF(xNew);
			if ( (fNew<fLo) || (fNew>fHi) )
				errMsg = sprintf( ...
				  "ERROR: Function is not monotonic: ( %g, %g, %g, %g, %g, %g ).", ...
				  xLo, xNew, xHi, fLo, fNew, fHi );
				error( errMsg );
			end
			if ( fNew < fTrgt )
				fLo = fNew;
				xLo = xNew;
			else
				fHi = fNew;
				xHi = xNew;
			end
			%
			if ( fLo < fMin )
				xBetter = xHi;
				fBetter = fHi;
			elseif ( fHi > fMax )
				xBetter = xLo;
				fBetter = fLo;
			elseif ( abs(fLo-fTrgt)<abs(fHi-fTrgt) )
				xBetter = xLo;
				fBetter = fLo;
			else
				xBetter = xHi;
				fBetter = fHi;
			end
			%
			subIter++;
			if ( abs(fBetter-fTrgt) <= fTol )
				break;
			elseif (subIter>=subIterLimit)
				break;
			end
		end % End sub loop.
		%
		numVals++;
		xVals(numVals) = xBetter;
		fVals(numVals) = fBetter;
	end % End main loop.
	numVals++;
	xVals(numVals) = x1;
	fVals(numVals) = f1;
	%
	dFVals = diff(fVals);
	minDF = min(dFVals);
	maxDF = max(dFVals);
	%
	if (flippedX)
		xVals_temp = xVals;
		xVals = xVals_temp(end:-1:1);
	end
	if (flippedF)
		fVals = -fVals;
	end
	msg_main( verbLev, thisFile, __LINE__, sprintf( ...
	  "Found %d points (out of %d) in %fs with interval range %9.2e to %9.2e.", ...
	  numVals, numValsRequested, time()-startTime, minDF, maxDF ) );
	retCode = RETCODE__SUCCESS;
	datOut.fVals = fVals;
return;
end

%!test
%!	commondefs;
%!	funchF = @(dummy)( (dummy.^0.1) + (dummy.^100.0) );
%!	x0 = 0.0;
%!	x1 = 1.0;
%!	numPts = 50;
%!	prm.verbLev = VERBLEV__MAIN;
%!	[ xVals, retCode ] = flinspace( x0, x1, numPts, funchF, prm );
%!	assert( RETCODE__SUCCESS == retCode );

%!test
%!	commondefs;
%!	funchF = @(dummy)( (dummy.^0.1) + (dummy.^100.0) );
%!	x0 = 1.0;
%!	x1 = 0.0;
%!	numPts = 50;
%!	prm.verbLev = VERBLEV__MAIN;
%!	[ xVals, retCode ] = flinspace( x0, x1, numPts, funchF, prm );
%!	assert( RETCODE__SUCCESS == retCode );

%!test
%!	commondefs;
%!	funchF = @(dummy)( -((dummy.^0.1) + (dummy.^100.0)) );
%!	x0 = 0.0;
%!	x1 = 1.0;
%!	numPts = 50;
%!	prm.verbLev = VERBLEV__MAIN;
%!	[ xVals, retCode ] = flinspace( x0, x1, numPts, funchF, prm );
%!	assert( RETCODE__SUCCESS == retCode );


%!test
%!	commondefs;
%!	funchF = @(dummy)( -((dummy.^0.1) + (dummy.^100.0)) );
%!	x0 = 1.0;
%!	x1 = 0.0;
%!	numPts = 50;
%!	prm.verbLev = VERBLEV__MAIN;
%!	[ xVals, retCode ] = flinspace( x0, x1, numPts, funchF, prm );
%!	assert( RETCODE__SUCCESS == retCode );
