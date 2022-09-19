function [ matJ, datOut ] = sja_scratch100( matV, matW, prm=[] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNZEPerRow = mygetfield( prm, "maxNZEPerRow", sizeK );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	%useVWeight = mygetfield( prm, "useVWeight", true );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNZEPerRow) );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
		%assert( isscalar(useVWeight) );
		%assert( islogical(useVWeight) );
	endif
	%
	% NZE = (identified) non-zero element.
	%
	matA = matV';
	for nf = 1 : sizeF
		msg( __FILE__, __LINE__, sprintf( "Analyzing row %d...", nf ) );
		vecB = matW(nf,:)';
		res0 = norm(vecB);
		%
		vecJ = zeros(sizeX,1);
		vecBeta = vecB - matA * vecJ;
		res = norm(vecBeta);
		nzeList = [];
		doRowLoop = true;
		while (doRowLoop)
			msg( __FILE__, __LINE__, sprintf( "  %3d:  %10.3f,  %10.3f.", max(size(nzeList)), res, norm(vecJ) ) );
			best_nx = 0;
			best_nzeList = nzeList;
			best_vecJ = vecJ;
			best_vecBeta = vecBeta;
			best_res = res;
			for trial_nx = 1 : sizeX
				if ( ismember( trial_nx, nzeList ) )
					continue;
				endif
				trial_nzeList = [ nzeList, trial_nx ];
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				if ( trial_res < best_res )
					best_nx = trial_nx;
					best_nzeList = trial_nzeList;
					best_vecJ = trial_vecJ;
					best_vecBeta = trial_vecBeta;
					best_res = trial_res;
				endif
			endfor
			nzeList = best_nzeList;
			vecJ = best_vecJ;
			vecBeta = best_vecBeta;
			res = best_res;
			clear trial_*;
			clear best_*;
			%
			if ( max(size(nzeList)) >= maxNZEPerRow )
				msg( __FILE__, __LINE__, sprintf( "Row %d: Reached maxNZEPerRow.", nf ) );
				doRowLoop = false; % Superfluous?
				break;
			elseif ( res < tol*res0 )
				msg( __FILE__, __LINE__, sprintf( "Row %d: Reached tol.", nf ) );
				doRowLoop = false; % Superfluous?
				break;
			endif
		endwhile
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
