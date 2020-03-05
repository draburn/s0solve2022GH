function [ vecV1, vecV2, matX, rvecD1, rvecD2, numD1Vals, numD2Vals ] = spanspace( ...
  vecX0, vecU1, vecU2, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealarray(vecU1,[sizeX,1]) );
	assert( isrealarray(vecU2,[sizeX,1]) );
	%
	s1Min = mygetfield( prm, "s1Min", -0.5 );
	s1Max = mygetfield( prm, "s1Max",  1.5 );
	s2Min = mygetfield( prm, "s2Min", -0.5 );
	s2Max = mygetfield( prm, "s2Max",  1.5 );
	assert( isrealscalar(s1Min) );
	assert( isrealscalar(s1Max) );
	assert( isrealscalar(s2Min) );
	assert( isrealscalar(s2Max) );
	assert( s1Min < s1Max );
	assert( s2Min < s2Max );
	%
	normU1Sq = vecU1'*vecU1;
	assert( 0.0 < normU1Sq );
	normU1 = sqrt(normU1Sq);
	vecV1 = vecU1 / normU1;
	vecT = vecU2 - (vecV1*(vecV1'*vecU2));
	normTSq = vecT'*vecT;
	assert( 0.0 < normTSq );
	normT = sqrt(normTSq);
	vecV2 = vecT / normT;
	%
	numD1Vals = mygetfield( prm, "numD1Vals", 51 );
	numD2Vals = mygetfield( prm, "numD2Vals", 53 );
	assert( isrealscalar(numD1Vals) );
	assert( isrealscalar(numD2Vals) );
	assert( 2 <= numD1Vals );
	assert( 2 <= numD2Vals );
	assert( round(numD1Vals) == numD1Vals );
	assert( round(numD2Vals) == numD2Vals );
	%
	d11 = vecV1'*vecU1;
	d12 = vecV1'*vecU2;
	d21 = vecV2'*vecU1; % Should be zero.
	d22 = vecV2'*vecU2;
	d1Lo = min([ 0.0, d11, d12 ]);
	d1Hi = max([ 0.0, d11, d12 ]);
	d2Lo = min([ 0.0, d21, d22 ]);
	d2Hi = max([ 0.0, d21, d22 ]);
	d1Mid = (d1Hi+d1Lo)/2.0;
	d1Var = (d1Hi-d1Lo)/2.0;
	d2Mid = (d2Hi+d2Lo)/2.0;
	d2Var = (d2Hi-d2Lo)/2.0;
	d1Min = d1Mid - d1Var + (s1Min*d1Var);
	d1Max = d1Mid + d1Var + (s1Max*d1Var);
	d2Min = d2Mid - d2Var + (s2Min*d2Var);
	d2Max = d2Mid + d2Var + (s2Max*d2Var);
	%
	%
	d1Vals = linspace( d1Min, d1Max, numD1Vals );
	d2Vals = linspace( d2Min, d2Max, numD2Vals );
	[ gridD1, gridD2 ] = meshgrid( d1Vals, d2Vals );
	rvecD1 = reshape( gridD1, 1, [] );
	rvecD2 = reshape( gridD2, 1, [] );
	matX = repmat(vecX0,[1,numD1Vals*numD2Vals]) + (vecV1*rvecD1) + (vecV2*rvecD2);
	%
return;
end


%!test
%!	tic();
%!	sizeX = 3;
%!	vecX0 = zeros(sizeX,1);
%!	vecU1 = randn(sizeX,1);
%!	vecU2 = randn(sizeX,1);
%!	[ vecV1, vecV2, matX, rvecD1, rvecD2, numD1Vals, numD2Vals ] = spanspace( vecX0, vecU1, vecU2 );
%!	numPts = size(matX,2);
%!	assert( numD1Vals*numD2Vals == numPts );
%!	for n=1:numPts
%!		rvecZ(n) = sqrt(min([ ...
%!		  sum((matX(:,n)-vecX0).^2), ...
%!		  sum((matX(:,n)-vecX0-vecU1).^2), ...
%!		  sum((matX(:,n)-vecX0-vecU2).^2) ]));
%!	end
%!	gridZ = reshape( rvecZ, numD2Vals, numD1Vals );
%!	gridX = reshape( rvecD1, numD2Vals, numD1Vals );
%!	gridY = reshape( rvecD2, numD2Vals, numD1Vals );
%!	contour( gridX, gridY, sqrt(gridZ), 31 );
%!	hold on;
%!	plot( [ 0.0, vecV1'*vecU1 ], [ 0.0, vecV2'*vecU1 ], 'ko-' );
%!	plot( [ 0.0, vecV1'*vecU2 ], [ 0.0, vecV2'*vecU2 ], 'kx-' );
%!	%axis equal;
%!	grid on;
%!	hold off;
%!	toc;
