function [ vecA, matV, fSS, vecGSS, matHSS, datOut ] = sshessfit1202( sizeX, numPts, matX, rvecF, matG, prm=[] )
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
	%%%msg( __FILE__, __LINE__, "USING GRAD IN V. THIS IS A BAD IDEA." );
	%%%matV = utorthdrop([ matDX(:,1:ptA-1), matDX(:,ptA+1:end), matG(:,ptA) ]);
	matV = utorthdrop([ matDX(:,1:ptA-1), matDX(:,ptA+1:end) ]);
	bigK = size(matV,2)
	matY = matV'*matDX;
	%%%matY = matV \ matDX;
	matGSS = matV'*matG;
	%
	hessfitPrm = [];
	%hessfitPrm.useCnstF = false;
	%hessfitPrm.useCnstG = false;
	%hessfitPrm.rTol = 1.0e-14;
	%hessfitPrm.maxIt = 2*sizeX;
	hessfitPrm = mygetfield( prm, "hessfitPrm", hessfitPrm );
	[ fSS, vecGSS, matHSS, hessfitDat ] = hessfit( bigK, numPts, matY, rvecF, matGSS, hessfitPrm );
	datOut.hessfitDat = hessfitDat;
return;
endfunction
