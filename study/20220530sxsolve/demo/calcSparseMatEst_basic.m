function [ matJEst, datOut ] = calcSparseMatEst_basic( matV, matW, prm = [] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARNING );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__LOW );
	maxNumElemPerRow = mygetfield( prm, "maxNumElemPerRow", floor(sizeK/2.0) );
	chiThresh = mygetfield( prm, "chiThresh", sqrt(eps) );
	useVWeight = mygetfield( prm, "useVWeight", true );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumElemPerRow) );
		assert( isrealscalar(chiThresh) );
		assert( 0.0 < chiThresh );
		assert( chiThresh < 1.0 );
		assert( isscalar(useVWeight) );
		assert( islogical(useVWeight) );
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
	matVVT = matV*(matV');
	matWVT = matW*(matV');
	if (useVWeight)
		matWeight = 1.0./(eps+1.0-matV2);
		rvecV2Avg = sum( matV2.*matWeight, 2 )';
		matVWeighted = matV.*matWeight;
	else
		rvecV2Avg = sum( matV2, 2 )';
	endif
	for m=1:sizeF
		cementedElementMsk = logical(zeros(1,sizeX));
		rvecJEst = zeros(1,sizeX);
		rvecRho = matW(m,:);
		matR = [];
		for l=1:maxNumElemPerRow
			if ( useVWeight )
				rvecRVAvg = sum( rvecRho.*matVWeighted, 2 )'; % Automatic broadcastig.
				rvecR2Avg = sum( (rvecRho.^2).*matWeight, 2 )'; % Automatic broadcasting.
			else
				rvecRVAvg = sum( rvecRho.*matV, 2 )'; % Automatic broadcastig.
				rvecR2Avg = sum( (rvecRho.^2), 2 )'; % Automatic broadcasting.
			endif
			rvecChi = (rvecRVAvg.^2)./( eps + rvecR2Avg.*rvecV2Avg );
			%
			% Find the element of rvecChi that is maximal among all uncemented elements
			%  and set the corresponding element of cementedElemMsk to true.
			[ chiMaxUncemented, indexOfMax ] = max( rvecChi(~cementedElementMsk) );
			if ( abs(chiMaxUncemented) < chiThresh )
				break;
			endif
			cementedElementMsk(( (1:sizeX)(~cementedElementMsk) )(indexOfMax) ) = true;
			%newElem = ( (1:sizeX)(~cementedElementMsk) )(indexOfMax);
			%cementedElementMsk( newElem ) = true;
			%
			%
			%%%rvecJEst(cementedElementMsk) = (matW(m,:)*(matV(cementedElementMsk,:)'))*inv(matV(cementedElementMsk,:)*(matV(cementedElementMsk,:)'));0)
			% The above rvecJEst is correct, but, let's optimize a bit...
			%
			matR = chol( matVVT(cementedElementMsk,cementedElementMsk) );
			% This cholinsert() seems to work, but, doesn't particularly help.
			%matR = cholinsert( matR, find(newElem==(1:sizeX)(cementedElementMsk)), matVVT(cementedElementMsk,newElem) );
			rvecJEst(cementedElementMsk) = (matWVT(m,cementedElementMsk)/matR)/(matR');
			%
			rvecRho = matW(m,:) - rvecJEst*matV;
		endfor
		matJEst(m,:) = rvecJEst;
	endfor
return;
endfunction
