function [ vecA, matV, fSS, vecGSS, matHSS, datOut ] = sshessfit1203( sizeX, numPts, matX, rvecF, matG, prm=[] )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	datOut = [];
	%
	[ foo, ptA ] = min(rvecF);
	vecA = matX(:,ptA);
	matDX = matX - vecA; % Autobroadcast.
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
	hessfitPrm = [];
	hessfitPrm = mygetfield( prm, "hessfitPrm", hessfitPrm );
	[ fSS, vecGSS, matHSS, hessfitDat ] = hessfit( bigK, numPts, matY, rvecF, matGSS, hessfitPrm );
	datOut.hessfitDat = hessfitDat;
return;
endfunction
