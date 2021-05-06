function xNext = pickNextPt1D_2pt( xVals, fVals, prm = [] )
  	% Should-be-precompiled...
	thisFile = "pickNextPoint1D_2pt.m";
	%
	% Check data types...
	assert(isrealarray(xVals));
	assert(2==size(size(xVals),2));
	assert(1==size(xVals,1));
	assert(isrealarray(fVals));
	assert(2==size(size(fVals),2));
	assert(1==size(fVals,1));
	assert(size(xVals,2)==size(fVals,2));
	numPts = size(xVals,2);
	%
	% Check requirements...
	assert(2==numPts);
	%
	% Check unsupported cases...
	assert(0==sum(0.0>=diff(xVals))); % Require xVals strictly increasing.
	avgF = sum(fVals)/numPts;
	avgFSq = sum(fVals.^2)/numPts;
	varSqF = avgFSq - (avgF^2);
	epsF = mygetfield(prm,"epsF",eps^0.5);
	assert(isrealarray(epsF,[1,1]));
	assert( 0.0 < epsF );
	assert(	varSqF > (epsF^2)*avgFSq ); % Require fVals not all same.
	assert(0==sum(0.0==fVals)); % Require fVals all non-zero.
	assert(0==sum(0.0>fVals(1)*fVals)); % Require all fVals have same sign.
	%
	% Do work.
	% Simple linear extrapolation.
	xNext = ( xVals(1)*fVals(2) - xVals(2)*fVals(1) ) / (fVals(2)-fVals(1));
	assert(isrealarray(xNext,[1,1]));
return;
end

%!test
%!	assert( 2.0 == pickNextPt1D_2pt( [0.0,1.0], [2.0,1.0] ) );
%!
%!	x0 = randn()
%!	x1 = x0 + abs(randn())
%!	f0 = randn()
%!	f1 = sign(f0)*abs(randn())
%!	xNext_mio = ( x0*f1 - x1*f0 ) / ( f1 - f0 )
%!	xNext = pickNextPt1D_2pt( [x0,x1], [f0,f1] )
%!	assert( abs(xNext-xNext_mio) < eps );
%!	
