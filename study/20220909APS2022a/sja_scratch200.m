function [ matJ, datOut ] = sja_scratch200( matV, matW, prm=[] )
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
		assert( maxNumNZEPerRow <= sizeK );
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
		while (doRowLoop)
			msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "  %3d:  %10.3f,  %10.3f.", max(size(nzeList)), norm(vecJ), res ) );
			%
			%msg( __FILE__, __LINE__, "cemented + 1..." );
			singlemenResVals = zeros(1,sizeX);
			for trial_nx = 1 : sizeX
				if ( ismember( trial_nx, nzeList ) )
					% We below assume these singlemenResVals(nx) are left at zero.
					continue;
				endif
				trial_nzeList = [ nzeList, trial_nx ];
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				if ( trial_res < tol*res0 )
					msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
					  "Row %d: (nx = %d) Short-circuit tol with NZE %d -> %d.", ...
					  nf, trial_nx, max(size(nzeList)), max(size(trial_nzeList)) ) );
					nzeList = trial_nzeList;
					doRowLoop = false;
					break;
				endif
				singlemenResVals(trial_nx) = trial_res;
			endfor
			if ( ~doRowLoop )
				clear trial_*;
				break;
			endif
			clear trial_*;
			%
			% Try short-circuit.
			if ( 1 )
				[ foo, trial_sorted ] = sort( singlemenResVals );
				trial_nzeList = trial_sorted(1:maxNumNZEPerRow);
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				if ( trial_res < tol*res0 )
					msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
					  "Row %d: Short-circuit tol with NZE %d -> %d.", nf, max(size(nzeList)), max(size(trial_nzeList)) ) );
					nzeList = trial_nzeList;
					doRowLoop = false;
					break;
				endif
			endif
			clear trial_*;
			%
			% That didn't work.
			% We'll probably need to cement one value.
			% But, which?
			% First, let's do a bit better than the singlemenResVals...
			%msg( __FILE__, __LINE__, "cemented + best + 1..." );
			[ foo, temp_sorted ] = sort( singlemenResVals );
			temp_nzeList = temp_sorted(1:maxNumNZEPerRow-1);
			% temp_nzeList is probably pretty good, so let's try it...
			% Let's loop over all elem and see if we can make it better.
			best_nx = 0;
			best_res = res0; % Any sufficiently large value is fine here.
			for trial_nx = 1 : sizeX
				if ( ismember( trial_nx, temp_nzeList ) )
					continue;
				endif
				trial_nzeList = [ temp_nzeList, trial_nx ];
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				if ( trial_res < tol*res0 )
					msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
					  "Row %d: (nx = %d) Short-circuit tol with NZE %d -> %d.", ...
					  nf, trial_nx, max(size(nzeList)), max(size(trial_nzeList)) ) );
					nzeList = trial_nzeList;
					doRowLoop = false;
					break;
				endif
				%
				if ( trial_res < best_res )
					best_nx = trial_nx;
					best_nzeList = trial_nzeList;
					best_vecJ = trial_vecJ;
					best_vecBeta = trial_vecBeta;
					best_res = trial_res;
				endif
			endfor
			if ( ~doRowLoop )
				clear trial_*;
				clear best_*;
				break;
			endif
			temp_nzeList = best_nzeList;
			clear trial_*;
			clear best_*;
			%
			%msg( __FILE__, __LINE__, "cemented + best - 1..." );
			% Great, now we need to decide which to cement.
			% The "best" is the worst one to drop.
			best_nx = 0;
			best_res = 0.0; % Any sufficiently small value is fine here.
			for trial_nx = 1 : sizeX
				if ( ismember( trial_nx, nzeList ) )
					continue;
				elseif ( ~ismember( trial_nx, temp_nzeList ) )
					continue;
				endif
				trial_nzeList = temp_nzeList( temp_nzeList ~= trial_nx );
				trial_vecJ = zeros(sizeX,1);
				trial_vecJ(trial_nzeList) = matA(:,trial_nzeList) \ vecB;
				trial_vecBeta = vecB - matA * trial_vecJ;
				trial_res = norm(trial_vecBeta);
				if ( trial_res < tol*res0 )
					% Note: this short-circuit should be unreachable in the current code.
					msgif( verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
					  "Row %d: (nx = %d) Short-circuit tol with NZE %d -> %d.", ...
					  nf, trial_nx, max(size(nzeList)), max(size(trial_nzeList)) ) );
					nzeList = trial_nzeList;
					doRowLoop = false;
					break;
				endif
				%
				if ( trial_res > best_res )
					best_nx = trial_nx;
					best_nzeList = trial_nzeList;
					best_vecJ = trial_vecJ;
					best_vecBeta = trial_vecBeta;
					best_res = trial_res;
				endif
			endfor
			if ( ~doRowLoop )
				clear trial_*;
				clear best_*;
				break;
			endif
			% Finally, cement this one.
			nzeList = [ nzeList, best_nx ];
			vecJ = zeros(sizeX,1);
			vecJ(nzeList) = matA(:,nzeList) \ vecB;
			vecBeta = vecB - matA * vecJ;
			res = norm(vecBeta);
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
			%msg( __FILE__, __LINE__, "Goodbye." );
			%return;
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
