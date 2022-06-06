function [ matJEst, datOut ] = calcSparseMatrixEstimate( matV, matW, prm = [] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__UNLIMITED );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__UNLIMITED );
	maxNumElemPerRow = mygetfield( prm, "maxNumElemPerRow", floor(sizeK/2.0) );
	chiThresh = mygetfield( prm, "chiThresh", sqrt(eps) );
	cEstRelThresh = mygetfield( prm, "cEstRelThresh", sqrt(eps) );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isposintscalar(sizeX) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumElemPerRow) );
		assert( isrealscalar(chiThresh) );
		assert( 0.0 < chiThresh );
		assert( chiThresh < 1.0 );
		assert( isrealscalar(cEstRelThresh) );
		assert( 0.0 < cEstRelThresh );
		assert( cEstRelThresh < 1.0 );
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
			
			%rvecWVAvg = sum( matW(m,:).*(matV.*matWeight), 2 )';
			%rvecW2Avg = sum( matW2(m,:).*matWeight, 2 )';
			%rvecChiW = (rvecWVAvg.^2)./( eps + rvecW2Avg.*rvecV2Avg );
			%sigmaChiW = sum(rvecChiW)
			
			%
			usedElemMsk = logical(zeros(1,sizeX));
			usedElemMsk(elemUsed) = true;
			selectorList = (1:sizeX)(~usedElemMsk);
			if (0)
				[ cEstMin, cEstOrderedList ] = sort( rvecCEst(~usedElemMsk), "descend" );
				[ chiMin, chiOrderedList ] = sort( rvecChi(~usedElemMsk), "descend" );
				% Or maybe add all that are above the other's other metric?
				if ( cEstMin < cEstRelThresh*max(rvecChi) )
					if ( chiMin < chiThresh )
						break;
					else
						elemUsed = [ elemUsed, selectorList(chiOrderedList(1)) ];
					endif
				else
					if ( chiMin < chiThresh )
						elemUsed = [ elemUsed, selectorList(cEstOrderedList(1)) ];
					elseif ( chiOrderedList(1) == cEstOrderedList(1) )
						elemUsed = [ elemUsed, selectorList(cEstOrderedList(1)) ];
					else
						elemUsed = [ elemUsed, selectorList(cEstOrderedList(1)), selectorList(chiOrderedList(1)) ];
					endif
				endif
			elseif (0)
				nBest = 0;
				trialNormBest = 0.0;
				for n = 1 : sizeX
				if ( sum(n==elemUsed)==0 )
					elemUsedTrial = [ elemUsed, n ];
					clear rvecJEstTrial;
					rvecJEstTrial = (matW(m,:)*(matV(elemUsedTrial,:)'))*inv(matV(elemUsedTrial,:)*(matV(elemUsedTrial,:)'));
					%trialNorm = norm(rvecJEstTrial,1);
					trialNorm = abs(rvecJEstTrial(end));
					if ( 0 == nBest || trialNorm < trialNormBest )
						nBest = n;
						trialNormBest = trialNorm;
					endif
				endif
				endfor
				elemUsed = [ elemUsed, nBest ];
			elseif (0)
				elemNotUsed = selectorList;
				clear rvecJOut;
				foo1 = (rvecRho(m,:)*(matV(elemNotUsed,:)')) ...
				 * inv( matV(elemNotUsed,:)*(matV(elemNotUsed,:)') + 1.0*sqrt(eps)*diag(diag(matV(elemNotUsed,:)*(matV(elemNotUsed,:)'))) );
				foo2 = (rvecRho(m,:)*(matV(elemNotUsed,:)')) ...
				 * inv( matV(elemNotUsed,:)*(matV(elemNotUsed,:)') + 2.0*sqrt(eps)*diag(diag(matV(elemNotUsed,:)*(matV(elemNotUsed,:)'))) );
				rvecJOut = 2.0*foo1 - foo2;
				%size(rvecJOut)
				%rvecJOut
				%selectorList
				[ foo, orderedListOfNotUsed ] = sort(abs(rvecJOut),"descend");
				newElem = selectorList(orderedListOfNotUsed(1));
				elemUsed = [ elemUsed, newElem ];
			else
				%[ foo, orderedList ] = sort( -rvecChi(~usedElemMsk) );
				[ foo, orderedList ] = sort( -abs(rvecCEst(~usedElemMsk)) );
				%[ foo, orderedList ] = sort( -abs(rvecCEst(~usedElemMsk).*rvecChi(~usedElemMsk)) );
				if ( abs(foo) < chiThresh )
					break;
				endif
				elemUsed = [ elemUsed, selectorList(orderedList(1)) ];
			endif
			%elemUsed
			
			rvecJEst(elemUsed) = (matW(m,:)*(matV(elemUsed,:)'))*inv(matV(elemUsed,:)*(matV(elemUsed,:)'));
			rvecRho = matW(m,:) - rvecJEst*matV;
			%norm(rvecRho,1)
		endfor
			
		useAugmentedList = true;
		if (useAugmentedList)
			if ( numel(elemUsed) > ceil(maxNumElemPerRow/2.0) )
				elemUsed = elemUsed( 1 : ceil(maxNumElemPerRow/2.0) );
				% Then, there's almost certainly no benefit to this aug list.
			endif
			%
			rvecJOut = zeros(1,sizeX);
			foo1 = (matW(m,:)*(matV(:,:)')) ...
			 * inv( matV(:,:)*(matV(:,:)') + 1.0*sqrt(eps)*diag(diag(matV(:,:)*(matV(:,:)'))) );
			foo2 = (matW(m,:)*(matV(:,:)')) ...
			 * inv( matV(:,:)*(matV(:,:)') + 2.0*sqrt(eps)*diag(diag(matV(:,:)*(matV(:,:)'))) );
			rvecJOut(:) = 2.0*foo1 - foo2;
			[ foo, altOrderedList ] = sort(abs(rvecJOut),"descend");
			%
			n = 0;
			while (numel(elemUsed)<maxNumElemPerRow)
				n++;
				if ( sum(altOrderedList(n)==elemUsed)==0 )
					elemUsed = [ elemUsed, altOrderedList(n) ];
				endif
			endwhile
			rvecJEst = zeros(1,sizeX);
			rvecJEst(elemUsed) = (matW(m,:)*(matV(elemUsed,:)'))*inv(matV(elemUsed,:)*(matV(elemUsed,:)'));
			%rvecJEst(elemUsed) = matW(m,:)/matV(elemUsed,:);
			%rvecJEst(elemUsed) = (matW(m,:)*(matV(elemUsed,:)'))*inv(matV(elemUsed,:)*(matV(elemUsed,:)')+1e-4*eye(maxNumElemPerRow));
			rvecRho = matW(m,:) - rvecJEst*matV;
			%sumsq(rvecRho)
			%rcond( matV(elemUsed,:)*(matV(elemUsed,:)') )
		endif
		elemUsed
			
		matJEst(m,:) = rvecJEst;
	endfor
return;
endfunction
