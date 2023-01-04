function [ matH, datOut ] = hessfit_simple_posdefy( matY, f0, vecGamma0, rvecF, matGamma, prm )
	datOut = [];
	%
	switch ( tolower(mygetfield( prm, "fitMethod", "" )) )
	case { "simple symm" }
		matA = (matY') \ (( matGamma - vecGamma0 )'); % Autobroadcast.
		matH = (matA'+matA)/2.0;
		return;
	case { "simple symm alt" }
		matA = matY'*(matGamma - vecGamma0); % Autobroadcast.
		matYTHY = (matA'+matA)/2.0;
		matH = matY' \ matYTHY / matY;
		matH = (matH'+matH)/2.0;
	case { "symm w diag" }
		%matH = (matA'+matA)/2.0 + (diag( vecADiag - diag(matA) ))/2.0;
		error( "Not implemented." );
	case { "", "posdefy loop" }
		vecYTGamma0 = matY' * vecGamma0;
		vecADiag = 2.0 * ( rvecF' - f0 - vecYTGamma0 );
		matA = matY'*matGamma - vecYTGamma0; % Autobroadcast.
		sz = size(matY,2);
		matYTHY = zeros(sz,sz)+inf;
		% Loops are slow in Octave, and eig should be unnecessary if matY is upper triangular. Perhaps optimize later.
		for n=1:sz
			for m=1:n-1
				if ( matA(m,n)*matA(n,m) <= 0.0 )
					matYTHY(n,m) = 0.0; % If it might be zero, make it zero.
				elseif ( abs(matA(m,n)) < abs(matA(n,m)) )
					matYTHY(n,m) = matA(m,n); % Take abs smaller.
				else
					matYTHY(n,m) = matA(n,m); % Take abs smaller.
				endif
				matYTHY(m,n) = matYTHY(n,m);
			endfor
			%matYTHY(n,n) = max([ matA(n,n), vecADiag(n) ]); % Take more positive value.
			if ( matA(n,n) > vecADiag(n) )
				matYTHY(n,n) = matA(n,n);
			else
				matYTHY(n,n) = vecADiag(n);
			endif
		endfor
		matH = matY' \ matYTHY / matY;
		matH = (matH'+matH)/2.0;
	otherwise
		error( "Invalid case." );
	endswitch
return
endfunction
