function matH = avech( vechL, sizeX=[], doChecks=false )
	if (isempty(sizeX))
		sizeX = floor(sqrt(2.0*size(vechL,1)));
	endif
	if (doChecks)
		assert( isposintscalar(sizeX) );	
		assert( isrealarray( vechL, [ (sizeX*(sizeX+1))/2, 1 ]) );
	endif
	%
	matL = zeros(sizeX,sizeX);
	m = 0;
	for n=1:sizeX
		matL(n:sizeX,n) = vechL(m+1:m+1+sizeX-n);
		m += 1+sizeX-n;
	endfor
	matH = matL' + matL;
return;
endfunction

%!test
%! sz = 5;
%! foo = randn(sz,sz);
%! matH0 = foo'*foo
%! matH0_mod = matH0;
%! for n = 1 : sz
%!	matH0_mod(n,n) /= 2.0;
%! endfor
%! vechL = vech(matH0_mod);
%! doChecks = true;
%! matH1 = avech( vechL, sz, doChecks )
%! res = sum(sum((matH1-matH0).^2))/sum(sum((matH0.^2)+(matH1.^2)))
%! assert( res < sqrt(eps) )
