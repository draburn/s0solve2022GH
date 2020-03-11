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
	marginSize = mygetfield( prm, "marginSize", 0.5 );
	%
	matV = myorth(matU);
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
	matVTU = matV'*matU;
	%
	for k=1:sizeK
		sA = min([ 0.0, matVTU(k,k) ]);
		sB = max([ 0.0, matVTU(k,k) ]);
		s0 = min([ 0.0, matVTU(k,:) ]);
		s1 = max([ 0.0, matVTU(k,:) ]);
		% Want values to span s0 ~ s1 plus some margin,
		% but, want the points sA and sB to be among the exact values.
		sMinIsh = s0 - (marginSize*(s1-s0));
		sMaxIsh = s1 + (marginSize*(s1-s0));
		numSVals = rvecN(k);
		vA = ceil(  numSVals*(sA-sMinIsh)/(sMaxIsh-sMinIsh) );
		vB = floor( numSVals*(sB-sMinIsh)/(sMaxIsh-sMinIsh) );
		if (vB<=vA)
			% We'll have to *force* things to be decent.
			vB = round( (vB+vA+1.0)/0.5 );
			vA = vB - 1;
		end
		assert( 1 <= vA );
		assert( vA < vB );
		assert( vB <= numSVals );
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
	% Crib from ndgrid to handle arbitrary case?
	numPts = prod( rvecN );
	matS = zeros( sizeK, numPts );
	switch (sizeK)
	case {1}
		matS(1,:) = dimDat(1).sVals;
	case {2}
		[ t1, t2 ] = ndgrid( dimDat(1).sVals, dimDat(2).sVals );
		matS(1,:) = reshape( t1, [1, numPts] );
		matS(2,:) = reshape( t2, [1, numPts] );
	case {3}
		[ matS1, matS2, matS3 ] = ndgrid( ...
		  dimDat(1).sVals, dimDat(2).sVals, dimDat(3).sVals );
		matS(1,:) = reshape( matS1, [1, numPts] );
		matS(2,:) = reshape( matS2, [1, numPts] );
		matS(3,:) = reshape( matS3, [1, numPts] );
	case {4}
		[ matS1, matS2, matS3, matS4 ] = ndgrid( ...
		  dimDat(1).sVals, dimDat(2).sVals, dimDat(3).sVals, dimDat(4).sVals );
		matS(1,:) = reshape( matS1, [1, numPts] );
		matS(2,:) = reshape( matS2, [1, numPts] );
		matS(3,:) = reshape( matS3, [1, numPts] );
		matS(4,:) = reshape( matS4, [1, numPts] );
	otherwise
		error(sprintf( "5+ dimensional case not supported." ));
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
