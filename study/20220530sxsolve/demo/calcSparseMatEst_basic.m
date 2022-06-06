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
	if (useVWeight)
		matWeight = 1.0./(eps+1.0-matV2);
		rvecV2Avg = sum( matV2.*matWeight, 2 )';
		matVWeighted = matV.*matWeight;
	else
		%matWeight = ones(sizeX,sizeK);
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
			% Find the element of rvecChi that is minimal among all uncemented elements.
			originalList = (1:sizeX)(~cementedElementMsk);
			[ chiMinUncemented, orderedList ] = sort( rvecChi(~cementedElementMsk), "descend" );
			if ( abs(chiMinUncemented) < chiThresh )
				break;
			endif
			newElem = originalList(orderedList(1));
			%
			if (cementedElementMsk(newElem))
				error( "INTERNAL ERROR: Attempting to add an element that has already been added. Please fix this bug!" );
			endif
			cementedElementMsk(newElem) = true;
			%
			%
			%%%rvecJEst(cementedElementMsk) = (matW(m,:)*(matV(cementedElementMsk,:)'))*inv(matV(cementedElementMsk,:)*(matV(cementedElementMsk,:)'));0)
			% The above rvecJEst is correct, but, let's optimize a bit...
			%
			matR = chol( matV(cementedElementMsk,:)*(matV(cementedElementMsk,:)') );
			% This cholinsert() seems to work, but, it feels icky.
			%%%matR = cholinsert( matR, find(newElem==(1:sizeX)(cementedElementMsk)), matV(cementedElementMsk,:)*(matV(newElem,:)') );
			rvecJEst(cementedElementMsk) = ((matW(m,:)*(matV(cementedElementMsk,:)'))/matR)/(matR');
			%
			%
			rvecRho = matW(m,:) - rvecJEst*matV;
		endfor
		matJEst(m,:) = rvecJEst;
	endfor
return;
endfunction
