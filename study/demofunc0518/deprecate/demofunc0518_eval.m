% 2D, support multi arg.

function matF = demofunc0518_eval( matX, funcPrm )
	numVals = size(matX,2);
	assert( isrealarray(matX,[2,numVals]) );
	assert( isposintscalar(funcPrm.numTerms1) || (0==funcPrm.numTerms1) );
	assert( isposintscalar(funcPrm.numTerms2) || (0==funcPrm.numTerms2) );
	%
	rvecR = sqrt(sum(matX.^2,1));
	rvecTheta = atan2(matX(2,:),matX(1,:));
	matF = zeros(2,numVals);
	%
	matF(1,:) = funcPrm.c10;
	if ( 1 <= funcPrm.numTerms1 )
		assert( isrealarray(funcPrm.c1,[funcPrm.numTerms1,4]) );
		for n=1:funcPrm.numTerms1
			matF(1,:) += funcPrm.c1(n,1) .* (rvecR.^funcPrm.c1(n,2)) ...
			  .* cos( (round(funcPrm.c1(n,3)) .* rvecTheta) + funcPrm.c1(n,4) );
		end
	end
	matF(2,:) = funcPrm.c20;
	if ( 1 <= funcPrm.numTerms2 )
		assert( isrealarray(funcPrm.c2,[funcPrm.numTerms2,4]) );
		for n=1:funcPrm.numTerms2
			matF(2,:) += funcPrm.c2(n,1) .* (rvecR.^funcPrm.c2(n,2)) ...
			  .* cos( (round(funcPrm.c2(n,3)) .* rvecTheta) + funcPrm.c2(n,4) );
		end
	end
end
