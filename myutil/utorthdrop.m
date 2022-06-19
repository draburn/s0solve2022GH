%  Function...
%   "ut" mean "upper-triangular".

function [ matV, rvecDrop ] = utorthdrop( matV, dropThresh = 1.0e-4 )
	sizeV = size(matV,2);
	rvecDrop = logical(zeros(1,sizeV));
	matV ./= ( eps + sqrt(sum(matV.^2,1)) );
	%
	vNorm = norm(matV(:,1));	
	if ( vNorm <= dropThresh )
		rvecDrop(1) = true;
		matV(:,1) = 0.0;
	endif
	%
	for k=2:sizeV
		matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
		vNorm = norm(matV(:,k));
		if ( vNorm <= dropThresh )
			rvecDrop(k) = true;
			matV(:,k) = 0.0;
		else
			matV(:,k) /= vNorm;
			matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
			vNorm = norm(matV(:,k));
			if ( vNorm <= dropThresh )
				rvecDrop(k) = true;
				matV(:,k) = 0.0;
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
%!	matV = utorthdrop_turbo(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();

%!test
%!	tic();
%!	sizeX = 10;
%!	sizeU = 50;
%!	matU = randn(sizeX,sizeU);
%!	[ matV, rvecDrop ] = utorthdrop_turbo(matU);
%!	sizeV = size(matV,2);
%!	assert( matV'*matV, eye(sizeV), (eps^0.75)*(sizeV^2) );
%!	toc();
