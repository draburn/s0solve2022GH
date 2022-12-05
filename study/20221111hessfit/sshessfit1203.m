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
	doOrtho = mygetfield( prm, "doOrtho", true );
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
	%
	inclGrad = mygetfield( prm, "inclGrad", false );
	if (inclGrad)
		error( "Not implemented!" );
		assert(doOrtho);
		vecU0 = matG(:,ptA);
		vecU = vecU0;
		vecU -= matV*(matV'*vecU);
		vecU -= matV*(matV'*vecU);
		dropThresh = sqrt(eps);
		if ( norm(vecU) > dropThresh * norm(vecU0) )
			vecU /= norm(vecU);
			matV = [ matV, vecU ];
			vecGSS = [vecGSS; 0.0];
			% Hmm.... maybe it should be included in solve??? Because H is sym?
			matHSS_old = matHSS;
			matHSS = zeros(bigK+1,bigK+1);
			matHSS(1:bigK,1:bigK) = matHSS_old;
			bigK++;
			%matHSS(bigK,bigK) = ???
		endif
	endif
return;
endfunction
