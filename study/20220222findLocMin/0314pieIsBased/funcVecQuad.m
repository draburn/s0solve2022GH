function vecF = funcVecQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0, prm = [] )
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	doTests = false;
	if (~isempty(prm))
		doTests = mygetfield( prm, "doTests", doTests );
	endif
	if (doTests)
		assert( isscalar(doTests) );
		assert( isbool(doTests) );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
	endif
	%
	vecDelta = vecX - vecX0;
	vecF = vecF0 + matJ0*vecDelta;
	for n=1:sizeF
		vecF(n) += 0.5*( vecDelta' * reshape( ary3Kappa0(n,:,:), [ sizeX, sizeX ] ) * vecDelta );
	endfor
endfunction

%!test
%!	setprngstates(0);
%!	sizeX = 4;
%!	sizeF = 4;
%!	vecX0 = randn(sizeX,1);
%!	vecF0 = randn(sizeF,1);
%!	matJ0 = randn(sizeF,sizeX);
%!	ary3Kappa0 = randn(sizeF,sizeX,sizeX);
%!	for n1=1:sizeX
%!	for n2=1:n1-1
%!		ary3Kappa0(:,n1,n2) += ary3Kappa0(:,n2,n1);
%!	endfor
%!	endfor
%!	for n1=1:sizeX
%!	for n2=n1+1:sizeX
%!		ary3Kappa0(:,n1,n2) = ary3Kappa0(:,n2,n1);
%!	endfor
%!	endfor
%!	%
%!	vecX = randn(sizeX,1);
%!	prm.doTests = true;
%!	funcVecQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0, prm )
