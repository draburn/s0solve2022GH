function matH = lambda2hess( vecLambda, datIn )
	matH = reshape( datIn.matDuplish * vecLambda, [datIn.sz,datIn.sz] );
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
