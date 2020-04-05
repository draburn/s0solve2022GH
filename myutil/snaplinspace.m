function  fVals = snaplinspace( fMin_approx, fMax_approx, numVals_approx, fSnapA, fSnapB )
	assert(isposintscalar(numVals_approx));
	assert(numVals_approx>=2);
	assert(isrealscalar(fMin_approx));
	assert(isrealscalar(fMax_approx));
	assert(isrealscalar(fSnapA));
	assert(isrealscalar(fSnapB));
	assert( fMin_approx <= fSnapA );
	assert( fSnapA < fSnapB );
	assert( fSnapB <= fMax_approx );
	%
	delta_approx = (fMax_approx-fMin_approx)/(numVals_approx-1.0);
	nA_approx = (fSnapA-fMin_approx)/delta_approx;
	nB_approx = (fSnapB-fMin_approx)/delta_approx;
	%
	%assert( nB_approx - nA_approx >= 1.0 );
	deltaN_approx = nB_approx - nA_approx;
	if ( deltaN_approx < 0.5 )
		f0 = (fSnapA+fSnapB)/2.0;
		fSnapA = f0 - (0.5*delta_approx);
		fSnapB = f0 + (0.5*delta_approx);
		deltaN_approx = 1;
	end
	%
	deltaN = round(deltaN_approx);
	nA = round(nA_approx);
	nB = nA + deltaN;
	assert( nB - nA >= 1 );
	delta = (fSnapB-fSnapA)/(nB-nA);
	numBelowA = ceil( (fSnapA-fMin_approx)/delta );
	numAboveB = ceil( (fMax_approx-fSnapB)/delta );
	%
	fMin = fSnapA - numBelowA*delta;
	fMax = fSnapB + numAboveB*delta;
	numVals = 1 + numBelowA + (nB-nA) + numAboveB;
	%
	fVals = linspace( fMin, fMax, numVals );
return;
end

%!test
%!	fMin_approx = 0.0;
%!	fMax_approx = 1.0;
%!	numVals_approx = 5;
%!	fSnapA = 0.1;
%!	fSnapB = 1.0/sqrt(2);
%!	fVals = snaplinspace( ...
%!	  fMin_approx, ...
%!	  fMax_approx, ...
%!	  numVals_approx, ...
%!	  fSnapA, ...
%!	  fSnapB );
%!	resA = min(abs(fVals-fSnapA));
%!	resB = min(abs(fVals-fSnapB));
%!	resX = min(abs(fVals-0.5));
%!	assert( resA < sqrt(eps) );
%!	assert( resB < sqrt(eps) );
%!	assert( resX > sqrt(eps) );
%!	assert( fVals(1) <= fMin_approx );
%!	assert( fVals(end) >= fMax_approx );

%!test
%!	fMin_approx = 0.0;
%!	fMax_approx = 1.0;
%!	numVals_approx = 5;
%!	fSnapA = 0.1;
%!	fSnapB = 1.0/sqrt(2);
%!	fSnapB *= 0.2;
%!	fVals = snaplinspace( ...
%!	  fMin_approx, ...
%!	  fMax_approx, ...
%!	  numVals_approx, ...
%!	  fSnapA, ...
%!	  fSnapB );
%!	resA = min(abs(fVals-fSnapA));
%!	resB = min(abs(fVals-fSnapB));
%!	assert( abs(resA-resB) <= eps );
%!	assert( fVals(1) <= fMin_approx );
%!	assert( fVals(end) >= fMax_approx );
