%  Function...
%    matV = myorth( matU, prm=[] )
%  Overview...
%    Generates an orthal matrix with the same span as a given matrix.
%    Similar to standard function ortho(), but, mine makes matV'*matU
%    be upper triangular.
%  Input values...
%   matU: The input matrix.
%   prm: Structure of parameters.
%  Output values...
%   matV: The output matrix.
%  See source code for information on prm.
function matV = myorth( matU, prm=[] )
	sizeX = size(matU,1);
	sizeK = size(matU,2);
	assert( isrealarray(matU,[sizeX,sizeK]) );
	assert( 1 <= sizeX );
	assert( 1 <= sizeK );
	assert( sizeK <= sizeX );
	%
	relTol = mygetfield( prm, "relTol", eps^1.5 );
	numPasses = mygetfield( prm, "numPasses", 2 );
	assert( isrealscalar(relTol) );
	assert( isrealscalar(numPasses) );
	assert( 0.0 < relTol );
	assert( 1 <= numPasses );
	%
	matV = matU;
	rvecUNormSq = sum( matU.^2, 1 );
	%
	k = 1;
	for n=1:numPasses
		f0 = matV(:,k)'*matV(:,k);
		if ( f0 <= relTol*rvecUNormSq(k) )
			error(sprintf( "Bad column vector at %d. (%g <= %g * %g).", ...
			  k, f0, relTol, rvecUNormSq(k) ));
		end
		matV(:,k) /= sqrt(f0);
	end
	%
	for k=2:sizeK
	for n=1:numPasses
		matV(:,k) -= matV(:,1:k-1) * ( matV(:,1:k-1)' * matV(:,k) );
		f0 = matV(:,k)'*matV(:,k);
		if ( f0 <= relTol*rvecUNormSq(k) )
			error(sprintf( "Bad column vector at %d during pass %d. (%g <= %g * %g).", ...
			  k, n, f0, relTol, rvecUNormSq(k) ));
		end
		matV(:,k) /= sqrt(f0);
	end
	end
return;
end

%!test
%!	tic();
%!	sizeX = 1000;
%!	sizeK = 100;
%!	matU = randn(sizeX,sizeK);
%!	matV = myorth(matU);
%!	assert( isrealarray(matV,[sizeX,sizeK]) );
%!	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
%!	toc();
