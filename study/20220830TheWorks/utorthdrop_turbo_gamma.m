%  Function...
%   "ut" mean "upper-triangular".

function [ matV, rvecDrop ] = utorthdrop_turbo_gamma( matV, dropThresh = 1.0e-4 )
	sizeV = size(matV,2);
	rvecDrop = logical(zeros(1,sizeV));
	matV ./= ( eps + sqrt(sum(matV.^2,1)) );
	%
	vNorm = norm(matV(:,1));	
	if ( vNorm <= dropThresh )
		rvecDrop(1) = true;
	endif
	%
	for k=2:sizeV
		matV(:,k) -= matV(:,~rvecDrop(1:k-1)) * ( matV(:,~rvecDrop(1:k-1))'*matV(:,k) );
		vNorm = norm(matV(:,k));
		if ( vNorm <= dropThresh )
			rvecDrop(k) = true;
		else
			matV(:,k) /= vNorm;
		endif		
	endfor
	for k=2:sizeV
	if (~rvecDrop(k))
		matV(:,k) -= matV(:,~rvecDrop(1:k-1)) * ( matV(:,~rvecDrop(1:k-1))'*matV(:,k) );
		vNorm = norm(matV(:,k));
		if ( vNorm <= dropThresh )
			rvecDrop(k) = true;
		else
			matV(:,k) /= vNorm;
		endif
	endif	
	endfor
	matV = matV(:,~rvecDrop);
return;
end

%!test
%!	tic();
%!	sizeX = 1000;
%!	sizeK = 100;
%!	matU = randn(sizeX,sizeK);
%!	matV = utorthdrop_turbo_gamma(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();

%!test
%!	tic();
%!	sizeX = 10;
%!	sizeU = 50;
%!	matU = randn(sizeX,sizeU);
%!	[ matV, rvecDrop ] = utorthdrop_turbo_gamma(matU);
%!	sizeV = size(matV,2);
%!	assert( matV'*matV, eye(sizeV), (eps^0.75)*(sizeV^2) );
%!	toc();
