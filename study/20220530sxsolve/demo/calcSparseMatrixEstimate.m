function [ matJEst, datOut ] = calcSparseMatrixEstimate( matV, matW, prm = [] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__UNLIMITED );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__UNLIMITED );
	maxNumElemPerRow = mygetfield( prm, "maxNumElem", floor(sizeK/2.0) );
	chiThresh = mygetfield( prm, "chiThresh", sqrt(eps) );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isposintscalar(sizeX) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumElemPerRow) );
		assert( isrealscalar(chiThresh) );
		assert( 0.0 < chiThresh );
		assert( chiThresh <= 1.0 );
	endif
	%
	matJEst = zeros(sizeF,sizeX); % Estimated sparse Jacobian.
	matCEst = zeros(sizeF,sizeX); % Like matJEst, but individually (not collectively) estimated coefficients.
	matChi = zeros(sizeF,sizeX); % Measures linear correlation of W with V.
	matRho = zeros(sizeF,sizeK); % Simply matW - matJEst * matV.
	datOut = [];
	%
	matV2 = matV.^2;
	matW2 = matW.^2;
	if (1)
		matWeight = 1.0./(eps+1.0-matV2);
	else
		matWeight = ones(sizeX,sizeK);
	endif
	for m=1:sizeF
		%
		elemUsed = [];
		rvecJEst = zeros(1,sizeX);
		rvecRho = matW(m,:);
		for l=1:maxNumElemPerRow
			rvecRVAvg = sum( rvecRho.*(matV.*matWeight), 2 )'; % Automatic broadcastig.
			rvecV2Avg = sum( matV2.*matWeight, 2 )';
			rvecR2Avg = sum( (rvecRho.^2).*matWeight, 2 )'; % Automatic broadcasting.
			rvecCEst = rvecRVAvg./rvecV2Avg;
			rvecChi = (rvecRVAvg.^2)./( eps + rvecR2Avg.*rvecV2Avg );
			%
			usedElemMsk = logical(zeros(1,sizeX));
			usedElemMsk(elemUsed) = true;
			selectorList = (1:sizeX)(~usedElemMsk);
			[ foo, orderedList ] = sort( -rvecChi(~usedElemMsk) );
			if ( abs(foo) < chiThresh )
				break;
			endif
			elemUsed = [ elemUsed, selectorList(orderedList(1)) ];
			rvecJEst(elemUsed) = (matW(m,:)*(matV(elemUsed,:)'))*inv(matV(elemUsed,:)*(matV(elemUsed,:)'));
			rvecRho = matW(m,:) - rvecJEst*matV;		
		endfor
		matJEst(m,:) = rvecJEst;
	endfor
return;
endfunction
