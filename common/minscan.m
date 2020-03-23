%  Function...
%    [ xBest, retCode, datOut ] = minscan( x0, x1, funchOmega, prm=[], datIn=[] )
%  Overview...
%    Searches for the minimum of funchOmega between x0 and x1.
%  Input values...
%    x0, x1: The bounding values in the argument for funchOmega.
%    funchOmega: A function handle from R^1 to R^1.
%      Negative values are probably allowed.
%    prm: Structure of parameters for the calculation.
%    datIn: Structure of additional input data.
%  Output values...
%    xBest: The value of x that best minimizes oemga.
%    retCode: A common return code, RETCODE__SUCCESS (0) on success.
%    datOut: Structure of additional output data.
%  See source code for more information on prm, datIn, and datOut.
%
%  This function is intended for use with postprocessing visualiation,
%   and may be inefficient and imprecise.
function [ xBest, retCode, datOut ] = minscan( x0, x1, funchOmega, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	commoninit;
	thisFile = "minscan";
	%
	assert( isrealscalar(x0) );
	assert( isrealscalar(x1) );
	assert( x0 < x1 );
	%
	numIterLimit = mygetfield(prm,"numIterLimit",5);
	numPtsPerIter = mygetfield(prm,"numPtsPerIter",11);
	funchOmegaSupportsMultiArg = mygetfield( prm, "funchOmegaSupportsMultiArg", false );
	assert( isposintscalar(numIterLimit) );
	assert( isposintscalar(numPtsPerIter) );
	assert( 3 <= numPtsPerIter );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	omega0 = funchOmega(x0);
	omega1 = funchOmega(x1);
	if (omega0<omega1)
		omegaBest = omega0;
		xBest = x0;
	else
		omegaBest = omega1;
		xBest = x1;
	end
	%
	xLo = x0;
	xHi = x1;
	for numIter=1:numIterLimit
		%
		rvecX = linspace( xLo, xHi, numPtsPerIter );
		if (funchOmegaSupportsMultiArg)
			rvecOmega = funchOmega(rvecX);
		else
			for n=1:numPtsPerIter
				rvecOmega(n) = funchOmega(rvecX(n));
			end
		end
		assert( isrealarray(rvecOmega,[1,numPtsPerIter]) );
		%
		[ omegaMin, indexOfMin ] = min(rvecOmega);
		if ( omegaMin < omegaBest )
			omegaBest = omegaMin;
			xBest = rvecX(indexOfMin);
		end
		if (indexOfMin==1)
			xHi = rvecX(2);
		elseif (indexOfMin==numPtsPerIter)
			xLo = rvecX(end-1);
		elseif (rvecOmega(indexOfMin+1)<rvecOmega(indexOfMin-1) )
			xLo = rvecX(indexOfMin);
			xHi = rvecX(indexOfMin+1);
		else
			xLo = rvecX(indexOfMin-1);
			xHi = rvecX(indexOfMin);
		end
	end
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_minscan
