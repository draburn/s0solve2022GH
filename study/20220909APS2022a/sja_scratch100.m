function [ matJ, datOut ] = sja_scratch100( matV, matW, prm=[] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNumNZEPerRow = mygetfield( prm, "maxNumNZEPerRow", sizeK-1 );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	%useVWeight = mygetfield( prm, "useVWeight", true );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumNZEPerRow) );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
		%assert( isscalar(useVWeight) );
		%assert( islogical(useVWeight) );
	endif
	%
	% NZE = (identified) non-zero element.
	%
	matA = matV';
	matJ = zeros(sizeF,sizeX);
	datOut = [];
	for nf = 1 : sizeF
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Analyzing row %d...", nf ) );
		vecB = matW(nf,:)';
		res0 = norm(vecB);
		%
		doRowLoop = true;
		nzeList = [];
		vecJ = zeros(sizeX,1);
		vecBeta = vecB - matA * vecJ;
		res = norm(vecBeta);
		%
		%vecG = matV*vecBeta;
		%g = norm(vecG);
		while (doRowLoop)
			msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  %3d:  %10.3f,  %10.3f.", max(size(nzeList)), norm(vecJ), res ) );
			best_nx = 0;
			best_nzeList = nzeList;
			best_vecJ = vecJ;
			best_vecBeta = vecBeta;
			best_res = res;
			%
			%best_vecG = vecG;
			%best_g = g;
			%msg( __FILE__, __LINE__, sprintf( "  nx %d yields res %10.3e and g %10.3e.", 0, res, g ) );
			%echo__vecJT = vecJ'
			%echo__vecBetaT = vecBeta'
			%echo__vecGT = vecG'
			%
			for trial_nx = 1 : sizeX
				if ( ismember( trial_nx, nzeList ) )
					continue;
				endif
				trial_nzeList = [ nzeList, trial_nx ];
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				%
				%trial_vecG = matV*trial_vecBeta;
				%trial_g = norm(trial_vecG); % Note: g = res because ||V*beta|| = ||beta||.
				%msg( __FILE__, __LINE__, sprintf( "  nx %d yields res %10.3e and g %10.3e.", trial_nx, trial_res, trial_g ) );
				%echo__trial_vecJT = trial_vecJ'
				%echo__trial_vecBetaT = trial_vecBeta'
				%echo__trial_vecGT = trial_vecG'
				%
				if ( trial_res < best_res )
					best_nx = trial_nx;
					best_nzeList = trial_nzeList;
					best_vecJ = trial_vecJ;
					best_vecBeta = trial_vecBeta;
					best_res = trial_res;
					%
					%best_vecG = trial_vecG;
					%best_g = trial_g;
				endif
			endfor
			nzeList = best_nzeList;
			vecJ = best_vecJ;
			vecBeta = best_vecBeta;
			res = best_res;
			clear trial_*;
			clear best_*;
			%
			if ( res < tol*res0 )
				msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
				  "Row %d: Reached tol with NZE %d.", nf, max(size(nzeList)) ) );
				doRowLoop = false; % Superfluous?
				break;
			elseif ( max(size(nzeList)) >= maxNumNZEPerRow )
				msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
				  "Row %d: Reached maxNumNZEPerRow with rel res %0.3e.", nf, res / res0 ) );
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
