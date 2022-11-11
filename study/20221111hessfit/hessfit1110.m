function [ f, vecG, matH, datOut ] = hessfit1110( matX, rvecF, matG, prm=[] )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	assert( 1 <= sizeX );
	assert( 1 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	rvecW1 = mygetfield( prm, "rvecW1", [] );
	if ( ~isempty(rvecW1) )
		assert( isrealarray(rvecW1,[1,numPts]) );
		[ sumS1, vecSumX1, matSumXXT1, vecC1, matC2 ] = __calcCnst_withW( matX, rvecF, matG, rvecW1 );
	else
		[ sumS1, vecSumX1, matSumXXT1, vecC1, matC2 ] = __calcCnst_sansW( matX, rvecF, matG );
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
	f = 0.0;
	vecG = zeros( sizeX, 1 );
	matH = zeros( sizeX, sizeX );
	datOut = [];
return;
endfunction


function [ sumW1, vecSumX1, matSumXXT1, vecC1, matC2 ] = __calcCnst_withW( matX, rvecF, matG, rvecW1 )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	assert( isrealarray(rvecW1,[1,numPts]) );
	%
	sumW1 = sum( rvecW1, 2 );
	%
	vecSumX1 = sum( rvecW1 .* matX, 2 ); % Autobroadcast rvecW1.
	%
	matSumXXT1 = zeros( sizeX, sizeX );
	for p=1:numPts
		matSumXXT1 += ( rvecW1(p) * matX(:,p) )' * matX(:,p);
	endfor
	%
	vecC1 = sum( rvecW1 .* matG, 2 ); % Autobroadcast rvecW1.
	%
	% Not sure about the most efficient way to calc matC2, but,
	%  this should be reasonable.
	matC2 = zeros( sizeX, sizeX );
	for p=1:numPts
		foo = ( rvecW1(p) .* matX(:,p) ) * ( matG(:,p)' ); % Autobroadcast rvecW1.
		matC2 += foo' + foo;
	endfor
return;
endfunction
function [ sumW1, vecSumX1, matSumXXT1, vecC1, matC2 ] = __calcCnst_sansW( matX, rvecF, matG )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	sumW1 = numPts;
	%
	vecSumX1 = sum( matX, 2 );
	%
	matSumXXT1 = zeros( sizeX, sizeX );
	for p=1:numPts
		matSumXXT1 += matX(:,p)' * matX(:,p);
	endfor
	%
	vecC1 = sum( matG, 2 );
	%
	% Not sure about the most efficient way to calc matC2,
	%  but, this should be reasonable.
	matC2 = zeros( sizeX, sizeX );
	for p=1:numPts
		foo = matX(:,p) * (matG(:,p)');
		matC2 += foo' + foo;
	endfor
return;
endfunction
