function [ matJ, datOut ] = sja_scratch000( matV, matW, prm=[] )
	mydefs;
	msg( __FILE__, __LINE__, "WARNING: This is discarded development code and probably does not work correctly." );
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNumLEPerRow = mygetfield( prm, "maxNumLEPerRow", sizeK );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	%useVWeight = mygetfield( prm, "useVWeight", true );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumLEPerRow) );
		%assert( isrealscalar(chiThresh) );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
		%assert( isscalar(useVWeight) );
		%assert( islogical(useVWeight) );
	endif
	%
	% LE = "large element"
	%    = non-zero element if matrix is truly sparse.
	%
	matVT = matV';
	%
	matD = matW;
	echo__sssqW = sum(sum(matW.^2))
	matL = zeros(sizeX,0); % List of LE for each row.
	doMainLoop = true;
	while (doMainLoop)
		[ foo, vecL ] = max( abs(matD*matVT), [], 2 ); % Tells us which elem to add to list for each row.
		% Above MIGHT lead to repetitions.
		% Worry about that later.
		matL = [ matL, vecL ]
		%echo__abs_matDVT = abs(matD*matVT)
		%echo__foo = foo
		%echo__rvecL = rvecL
		%echo__matL = matL
		matJ = zeros(sizeF,sizeX);
		for m = 1 : sizeF
			matJ(m,matL(m,:)) = matD(m,:) / matV(matL(m,:),:);
		endfor
		matD = matW - (matJ * matV);
		%echo__matJ = matJ
		%echo__matD = matD
		%echo__size_matL_2 = size(matL,2)
		%echo__maxNumLEPerRow = maxNumLEPerRow
		echo__sssqD = sum(sum(matD.^2))
		%echo__sssqW = sum(sum(matW.^2))
		%
		if ( size(matL,2) >= maxNumLEPerRow )
			msg( __FILE__, __LINE__, "Reached maxNumLEPerRow." );
			doMainLoop = false; % Superfluous?
			break;
		elseif ( sum(sum(matD.^2)) < tol*sum(sum(matW)) )
			msg( __FILE__, __LINE__, "Reached tol." );
			doMainLoop = false; % Superfluous?
			break;
		endif
	endwhile
	datOut = [];
return;
endfunction
