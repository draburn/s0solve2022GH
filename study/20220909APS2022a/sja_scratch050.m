function [ matJ, datOut ] = sja_scratch050( matV, matW, prm=[] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNumNZEPerRow = mygetfield( prm, "maxNumNZEPerRow", sizeK-1 );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumNZEPerRow) );
		assert( maxNumNZEPerRow <= sizeK );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
	endif
	%
	% NZE = (identified) non-zero element.
	matA = matV';
	matJ = zeros(sizeF,sizeX);
	datOut = [];
	for nf = 1 : sizeF
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Analyzing row %d...", nf ) );
		vecB = matW(nf,:)';
		res0 = norm(vecB);
		if ( 0.0 == res0 )
			matJ(nf,:) = vecJ;
			continue;
		endif
		%
		doRowLoop = true;
		minListSize = 0; % row loop iterations.
		nzeList = [];
		vecJ = zeros(sizeX,1);
		vecBeta = vecB - matA * vecJ;
		res = norm(vecBeta);
		prev_nzeList = [];
		while (doRowLoop)
			msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  %3d:  %10.3f,  %10.3f.", max(size(nzeList)), norm(vecJ), res ) );
			%
			% Look at [ pseudo-cemented, * ].
			nzeMask = logical(zeros(1,sizeX));
			nzeMask(nzeList) = true;
			singlemenResVals = zeros(1,sizeX); % Code below counts on this suffic small for the nzeList elements.
			for trial_nx = (1:sizeX)(~nzeMask)
				singlemenResVals(trial_nx) = norm( vecB - matA(:,[nzeList,trial_nx]) * ( matA(:,[nzeList,trial_nx]) \ vecB ) );
			endfor
			%
			% Check short-circuit for [ pseudo-cemented, several best per signlemen ].
			[ foo, temp_sorted ] = sort( singlemenResVals ); % Assumes singlemenResVals(nzeList) is suffic small (0.0).
			temp_nzeList = temp_sorted(1:maxNumNZEPerRow);
			if ( verbLev >= VERBLEV__INFO )
				msg( __FILE__, __LINE__, "After singlemen analysis..." );
				temp_nzeList
			endif
			temp_vecJ = zeros(sizeX,1);
			temp_vecJ(temp_nzeList) = matA(:,temp_nzeList) \ vecB;
			temp_res = norm( vecB - matA * temp_vecJ );
			if ( temp_res < tol*res0 )
				msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( "Row %d: Reached tol with minListSize %d.", nf, minListSize ) );
				nzeList = temp_nzeList;
				doRowLoop = false; % Superfluous?
				break;
			elseif ( minListSize+1 >= maxNumNZEPerRow )
				msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( "Row %d: Reached maxNumNZEPerRow with rel res %0.3e.", nf, res / res0 ) );
				nzeList = temp_nzeList;
				doRowLoop = false; % Superfluous?
				break;
			elseif ( ~isempty(prev_nzeList) && 0 == sum( temp_nzeList ~= prev_nzeList ) )
				% We could actually perform this check before the multiemen loop and fetch the result from the previous iter,
				% but, I doubt that would be worth the effort.
				msgif( verbLev >= VERBLEV__FLAG, __FILE__, __LINE__, sprintf( ...
				  "Row %d: Reached repetition with rel res %0.3e and minListSize %d.", nf, best_res / res0, minListSize ) );
				nzeList = temp_nzeList;
				doRowLoop = false; % Superfluous?
				break;
			endif
			prev_nzeList = temp_nzeList;
			if ( verbLev >= VERBLEV__INFO )
				msg( __FILE__, __LINE__, "After setting...." );
				nzeList
			endif
			%
			minListSize++;
			nzeList = temp_sorted(1:minListSize);
		endwhile
		vecJ = zeros(sizeX,1);
		vecJ(nzeList) = matA(:,nzeList) \ vecB;
		vecBeta = vecB - matA * vecJ;
		res = norm(vecBeta);
		%
		matJ(nf,:) = vecJ;
		clear vecB;
		clear res0;
		clear vecJ;
		clear vecBeta;
		clear res;
		clear nzeList;
		clear doRowLoop;
	endfor
	%
	datOut = [];
	return;
return;
endfunction
