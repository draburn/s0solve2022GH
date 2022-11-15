function [ vecX0, matV, f, vecGss, matHss, datOut ] = sshessfit1115_nonorth( sizeX, numPts, matX, rvecF, matG, pt0, prm=[] )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	assert( isposintscalar(pt0) );
	assert( pt0 <= numPts );
	assert( sizeX > numPts );
	%
	vecX0 = matX(:,pt0);
	matV = [ matX(:,1:pt0-1), matX(:,pt0+1:end) ] - vecX0;
	bigK = numPts-1;
	matI = eye(bigK);
	matY = [ matI(:,1:pt0-1), zeros(bigK,1), matI(:,pt0:end) ];
	% Leaving the zero in place is eaiser, though perhaps less elegant.
	matGss = matV'*matG; % ... Yeah?
	%
	hessfitPrm = mygetfield( prm, "hessfitPrm", [] );
	[ f, vecGss, matHss, hessfitDat ] = hessfit( bigK, numPts, matY, rvecF, matGss, hessfitPrm );
	datOut = [];
	datOut.hessfitDat = hessfitDat;
return;
endfunction
