function [ vecA, matV, f, vecG, matH, datOut ] = sshessfit( matX, rvecF, matG, prm=[] )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( 2 <= sizeX );
	assert( 2 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	ptA = mygetfield( prm, "ptA", 1 );
	assert( isposintscalar(ptA) );
	assert( ptA <= numPts );
	vecA = matX(:,ptA);
	matDX = matX - vecA; % Autobroadcast.
	%
	if ( numPts >= sizeX+1 )
		msg( __FILE__, __LINE__, "Note: numPts >= sizeX+1. Calling hessfit on full space." );
		matV = eye(sizeX,sizeX);
		[ f, vecG, matH, datOut ] = hessfit( matDX, rvecF, matG, prm );
		return;
	endif
	%
	doOrtho = mygetfield( prm, "doOrtho", false );
	if (doOrtho)
		matV = utorthdrop([ matDX(:,1:ptA-1), matDX(:,ptA+1:end) ]);
		bigK = size(matV,2);
		matY = matV'*matDX;
	else
		matV = [ matDX(:,1:ptA-1), matDX(:,ptA+1:end) ];
		bigK = size(matV,2);
		%matY = matV \ matDX;
		matI = eye(bigK,numPts-1);
		matY = [ matI(:,1:ptA-1), zeros(bigK,1), matI(:,ptA:end) ]; % Right?
	endif
	matGSS = matV'*matG;
	%
	[ f, vecG, matH, datOut ] = hessfit( matY, rvecF, matGSS, prm );
return;
endfunction
