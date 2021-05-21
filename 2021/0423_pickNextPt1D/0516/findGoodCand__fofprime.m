function [ xCand, meritCand ] = findGoodCand__fofprime( ...
  xVals_input, gVals_input, prm = [] );
  	% Should-be-precompiled...
	commondefs;
	thisFile = "findGoodCand__fofprime";
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	xCand = [];
	meritCand = -1.0;
	%return;
	%
	numPts = size(xVals_input,2);
	% We'll reverse all of these to test for the right side.
	xVals = xVals_input;
	gVals = gVals_input;
	dxVals = diff(xVals);
	cxVals = cent(xVals);
	dgVals = diff(gVals);
	cgVals = cent(gVals);
	chVals = cgVals .* dxVals ./ dgVals;
	%
	numPtsForFit = mygetfield( prm, "numPtsForFit", 2 );
	numPtsToCheck = mygetfield( prm, "numPtsToCheck", 4 );
	assert( numPtsForFit >= 2 );
	assert( numPtsToCheck >= 3 );
	assert( numPtsToCheck >= numPtsForFit );
	%
	for n = numPtsToCheck : numPts-1
	if ( gVals(n) < gVals(n-1) && gVals(n) < gVals(n+1) )
		meritCand = -1.0;
		findGoodCand__fofprime__sub;
		thisFile = "findGoodCand__fofprime";
		if ( 0.0 < meritCand )
			return;
		end
	end
	end
	%
	xCand = [];
	meritCand = -1.0;
return;
end
