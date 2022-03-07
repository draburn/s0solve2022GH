%  Function...

function [ matV, rvecDrop ] = utorthdrop( matU, tol=0.9999, prm=[] )
	sizeX = size(matU,1);
	sizeU = size(matU,2);
	assert( isrealarray(matU,[sizeX,sizeU]) );
	%
	numPasses = mygetfield( prm, "numPasses", 2 );
	assert( isrealscalar(tol) );
	assert( isrealscalar(numPasses) );
	assert( 0.0 < tol );
	assert( tol < 1.0 );
	assert( 1 <= numPasses );
	
	
	rvecDrop = logical(zeros(1,sizeU));
	matV = matU;
	sizeV = sizeU;
	for n=1:numPasses
		for k=1:sizeU
		if (~rvecDrop(k))
			vSqBefore = sumsq(matV(:,k));
			if ( k > 1 )
				% May involve useless calc for zero vectors,
				%  but, implementation is easier this way.
				matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)'*matV(:,k) );
			endif
			vSqAfter = sumsq(matV(:,k));
			if ( vSqAfter <= (1.0-tol)*vSqBefore )
				rvecDrop(k) = true;
				matV(:,k) = 0.0;
			else
				assert( vSqAfter > 0.0 );
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
%!	matV = utorthdrop(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();
