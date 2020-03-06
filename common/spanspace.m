function [ matV, matX, matS ] = spanspace( matU, rvecN, prm=[] )
	sizeX = size(matU,1);
	sizeK = size(matU,2);
	assert( isrealarray(matU,[sizeX,sizeK]) );
	assert( isrealarray(rvecN,[1,sizeK]) );
	assert( sizeK >= 1 );
	assert( sizeK <= sizeX );
	assert( rvecN(:) >= 1 );
	assert( round(rvecN(:)), rvecN(:), -(eps^0.75) );
	%
	matV = myorth(matU);
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
	matVTU = matV'*matU;
	%
	for k=1:sizeK
		sA = min([ 0.0, matVTU(k,1:k) ]);
		sB = max([ 0.0, matVTU(k,1:k) ]);
		s0 = min([ 0.0, matVTU(k,:) ]);
		s1 = max([ 0.0, matVTU(k,:) ]);
		% Want values to span s0 ~ s1 plus some margin,
		% but, want the points sA and sB to be among the exact values.
		sMinIsh = s0 - (0.2*(s1-s0));
		sMaxIsh = s1 + (0.2*(s1-s0));
		numSVals = rvecN(k);
		vA = ceil(  numSVals*(sA-sMinIsh)/(sMaxIsh-sMinIsh) );
		vB = floor( numSVals*(sB-sMinIsh)/(sMaxIsh-sMinIsh) );
		assert( vB > vA );
		dsdv = (sB-sA)/(vB-vA);
		sMin = sA - (dsdv*(vA-1));
		sMax = sB + (dsdv*(numSVals-vB));
		sVals = linspace( sMin, sMax, numSVals );
		assert( sVals(vA), sA, (eps^0.75)*(sMax-sMin) );
		assert( sVals(vB), sB, (eps^0.75)*(sMax-sMin) );
		%
		dimDat(k).numSVals = numSVals;
		dimDat(k).sVals = sVals;
	end
	%
	% We could probably instead use ndgrid... If I understood how it worked.
	numPts = prod( rvecN );
	matS = zeros( sizeK, numPts );
	switch (sizeK)
	case {1}
		matS(1,:) = dimDat(1).sVals;
	case {2}
		[ t1, t2 ] = meshgrid( dimDat(1).sVals, dimDat(2).sVals );
		matS(1,:) = reshape( t1, [1, numPts] );
		matS(2,:) = reshape( t2, [1, numPts] );
	case {3}
		error(sprintf( "3 dimensional case not supported." ));
	otherwise
		error(sprintf( "4+ dimensional case not supported." ));
	end
	%
	matX = matV * matS;
return;
end

%!test
%!	sizeX = 20;
%!	sizeK = 2;
%!	matU = randn(sizeX,sizeK);
%!	rvecN = 5*ones(1,sizeK);
%!	[ matV, matX, matS ] = spanspace( matU, rvecN );
