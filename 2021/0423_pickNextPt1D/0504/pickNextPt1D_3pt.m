function xNext = pickNextPt1D_3pt( xVals, fVals, prm = [] )
  	% Should-be-precompiled...
	thisFile = "pickNextPoint1D_3pt.m";
	%
	% Check data types...
	numPts = size(xVals,2);
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	%
	% Check requirements...
	assert( 3 == numPts );
	%
	% Check unsupported cases...
	xValsAreStrictlyIncreasing =( 0==sum(0.0>=diff(xVals)) );
	assert(xValsAreStrictlyIncreasing);
	%
	avgF = sum(fVals)/numPts;
	avgFSq = sum(fVals.^2)/numPts;
	varSqF = avgFSq - (avgF^2);
	epsF = mygetfield(prm,"epsF",eps^0.5);
	assert( isrealarray(epsF,[1,1]) );
	assert( 0.0 < epsF );
	fValRangeIsSufficient =( varSqF > (epsF^2)*avgFSq );
	assert( fValRangeIsSufficient );
	%
	signF = sign(fVals(1));
	gVals = signF * fVals;
	fValsAllHaveSameSign =( 0 == sum(0.0>=gVals) );
	assert( fValsAllHaveSameSign );
	%
	% Do work.
	error( "Not implemented." );
	%xNext = 2.0*xVals(3) - xVals(1);
	%assert( isrealarray(xNext,[1,1]) );
return;
end

%!test
%!	x1 = randn()
%!	x2 = x1 + abs(randn())
%!	x3 = x2 + abs(randn())
%!	f1 = randn()
%!	f2 = sign(f1)*abs(randn())
%!	f3 = sign(f1)*abs(randn())
%!	xNext = pickNextPt1D_3pt( [x1,x2,x3], [f1,f2,f3] )
%!	
