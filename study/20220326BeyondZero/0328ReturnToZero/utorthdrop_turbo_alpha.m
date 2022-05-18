%  Function...
%   "ut" mean "upper-triangular".

function [ matV, rvecDrop ] = utorthdrop_turbo_alpha( matV, tol=0.9999 )
	sizeV = size(matV,2);
	rvecDrop = logical(zeros(1,sizeV));
	for n=1:2
		for k=1:sizeV
		if (~rvecDrop(k))
			vSqBefore = sumsq(matV(:,k));
			if ( k > 1 )
				matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
			endif
			vSqAfter = sumsq(matV(:,k));
			if ( vSqAfter <= (1.0-tol)*vSqBefore )
				rvecDrop(k) = true;
				matV(:,k) = 0.0;
			else
				matV(:,k) /= sqrt(vSqAfter);
			endif
		endif
		endfor
	endfor
	matV = matV(:,~rvecDrop);
return;
end

%!test
%!	tic();
%!	sizeX = 1000;
%!	sizeK = 100;
%!	matU = randn(sizeX,sizeK);
%!	matV = utorthdrop_turbo_alpha(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();
