function [ vecLambda, datOut ] = hess2lambda( matHess )
	vecLambda = vech( matHess - diag(diag(matHess)/2.0) );
	if ( 1 == nargout )
		return;
	endif
	sz = size(matHess,1);
	datOut.sz = sz;
	datOut.matDuplish = sparse(duplication_matrix(sz));
	% We need to halve the coeff corresponding to the diagona of the hessian.
	for n = 1 : sz
		m = (n-1)*sz + n;
		datOut.matDuplish( m, m - (n*(n-1))/2 ) = 2.0;
	endfor
return;
endfunction

% matH = matL' + matL
% vecLambda = vech(matL)
%  where "vech" represents "half-vectorization".
% matH = reshape( matDuplish*vecLambda, [sz,sz] )
%  where matDuplish is like a "duplication matrix" (but not exactly).

%!test
%!	sz = 5;
%!	foo = randn(sz,sz);
%!	matH0 = foo'*foo
%!	[ vecLambda, dat ] = hess2lambda( matH0 );
%!	matH1 = lambda2hess( vecLambda, dat )
%!	matR = matH1 - matH0
%!	res = sum(sum(matR.^2))/sum(sum((matH0.^2)+(matH1.^2)))
%!	assert( res < sqrt(eps) )
