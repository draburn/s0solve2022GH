%  Function...
%   "ut" mean "upper-triangular".

function [ matV ] = utorth_turbo( matV, tol=0.9999 )
	matV(:,1) /= norm(matV(:,1));
	for k=2:size(matV,2)
		matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
		matV(:,k) /= norm(matV(:,k));
	endfor
	for k=2:size(matV,2)
		matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
		matV(:,k) /= norm(matV(:,k));
	endfor
return;
end

%!test
%!	tic();
%!	sizeX = 1000;
%!	sizeK = 100;
%!	matU = randn(sizeX,sizeK);
%!	matV = utorth_turbo(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();
