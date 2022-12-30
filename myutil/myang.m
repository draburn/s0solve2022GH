function ang = myang( vecA, vecB );
	sz = size(vecA,1);
	assert( isrealarray(vecA,[sz,1]) );
	assert( isrealarray(vecB,[sz,1]) );
	nA = norm(vecA);
	nB = norm(vecB);
	if ( 0.0 == nA || 0.0 == nB )
		ang = 0.0;
		return;
	endif
	ang = (vecA'*vecB)/(nA*nB);
return;
end
