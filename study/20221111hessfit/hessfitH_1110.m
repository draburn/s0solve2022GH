function [ matH, datOut ] = hessfitH_1110( sizeX, numPts, vecG, matX, matG, prm=[] )
	error( "Yo! This is a work-in-progress!" );
	
	assert( 1 <= sizeX );
	assert( 1 <= numPts );
	assert( isrealarray(vecG,[sizeX,1]) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	rvecW1 = mygetfield( prm, "rvecW1", [] );
	if ( ~isempty(rvecW1) )
		assert( isrealarray(rvecW1,[1,numPts]) );
		[ sumW1, vecSumX1, matSumXXT1, matC2 ] = __calcCnst_withW1( sizeX, numPts, matX, matG, rvecW1 );
	else
		[ sumW1, vecSumX1, matSumXXT1, matC2 ] = __calcCnst_sansW1( sizeX, numPts, matX, matG );
	endif
	clear rvecW1;
	%
	epsHRegu = mygetfield( prm, "epsHRegu", [] );
	rvecW2 = mygetfield( prm, "rvecW2", [] );
	if ( ~isempty(rvecW2) )
		if ( ~isempty(epsHReg) )
			error( "prm.epsHRegu and prm.rvecW2 were both specified." );
		else
			assert( isrealarray(rvecW2,[1,numPts]) );
			epsHReg = 2.0 * sum( rvecW2, 2 );
		endif
	else
		if ( isempty(epsHRegu) )
			epsHRegu = 0.0; % This is maybe safe?
		endif
	endif
	clear rvecW2;
	%
	%
	%
	matH_initialGuess = mygetfield( prm, "matH_initialGuess", zeros(sizeX,sizeX) );
	assert( isrealarray(matH_initialGuess,[sizeX,sizeX]) );
	assert( issymmetric(matH_initialGuess) );
	%
	%
	%
	vecL_initialGuess = vech(matH_initialGuess);
	
	
	msg( __FILE__, __LINE__, "Begin test code..." );
	vecResL = __funcResL( sizeX, numPts, matH_initialGuess, vecG, sumW1, vecSumX1, matSumXXT1, matC2, epsHRegu );
	size__vecReL = size(vecResL)
	echo__normResL = norm(vecResL)
	
	
	matH = zeros( sizeX, sizeX );
	datOut = [];
return;
endfunction


function [ sumW1, vecSumX1, matSumXXT1, matC2 ] = __calcCnst_withW1( sizeX, numPts, matX, matG, rvecW1 )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	assert( isrealarray(rvecW1,[1,numPts]) );
	%
	sumW1 = sum( rvecW1, 2 );
	vecSumX1 = sum( rvecW1 .* matX, 2 ); % Autobroadcast rvecW1.
	matSumXXT1 = zeros( sizeX, sizeX );
	for p=1:numPts
		matSumXXT1 += ( rvecW1(p) * matX(:,p) )' * matX(:,p);
	endfor
	%
	matC2 = zeros( sizeX, sizeX );
	for p=1:numPts
		foo = ( rvecW1(p) .* matX(:,p) ) * ( matG(:,p)' ); % Autobroadcast rvecW1.
		matC2 += foo' + foo;
	endfor
	% Not sure about the most efficient way to calc matC2, but, this should be reasonable.
return;
endfunction
function [ sumW1, vecSumX1, matSumXXT1, matC2 ] = __calcCnst_sansW1( sizeX, numPts, matX, matG )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	sumW1 = numPts;
	vecSumX1 = sum( matX, 2 );
	matSumXXT1 = zeros( sizeX, sizeX );
	for p=1:numPts
		matSumXXT1 += matX(:,p)' * matX(:,p);
	endfor
	%
	matC2 = zeros( sizeX, sizeX );
	for p=1:numPts
		foo = matX(:,p) * (matG(:,p)');
		matC2 += foo' + foo;
	endfor
	% Not sure about the most efficient way to calc matC2, but, this should be reasonable.
return;
endfunction


function vecResL = __funcResL( sizeX, numPts, matH, vecG, sumW1, vecSumX1, matSumXXT1, matC2, epsHRegu )
	foo = vecG * (vecSumX1') + matH * matSumXXT1;
	vecResL = vech( foo' + foo - matC2 + (epsHRegu*matH) );
return;
endfunction
