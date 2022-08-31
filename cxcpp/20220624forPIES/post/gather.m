clear;

filepathBase = "./00000000/runs/jacobian__82445x31/";
vecX0 = load( [filepathBase "2374802/jevalDriverOutput/x000000.txt"] );
vecF0 = load( [filepathBase "2374802/jevalDriverOutput/f000000.txt"] );
assert( isreal(vecX0) );
assert( isreal(vecF0) );

sizeX = size(vecX0,1);
sizeF = size(vecF0,1);
assert( size(vecX0,2) == 1 );
assert( size(vecF0,2) == 1 );

sizeK = sizeX;
rvecEpsX = zeros(1,sizeK);
matY = zeros(sizeX,sizeK);
matJ = zeros(sizeF,sizeK);
vecFSumSq = vecF0.^2;

for m=0:7
	switch (m)
	case 0
		basePath = [filepathBase "2374802/jevalDriverOutput"];
	case 1
		basePath = [filepathBase "2374803/jevalDriverOutput"];
	case 2
		basePath = [filepathBase "2374804/jevalDriverOutput"];
	case 3
		basePath = [filepathBase "2381797/jevalDriverOutput"];
	case 4
		basePath = [filepathBase "2381798/jevalDriverOutput"];
	case 5
		basePath = [filepathBase "2381799/jevalDriverOutput"];
	case 6
		basePath = [filepathBase "2383275/jevalDriverOutput"];
	case 7
		basePath = [filepathBase "2383276/jevalDriverOutput"];
	otherwise
		error( "Invalid case." );
	end
	for n=1:100
		k = 100*m + n;
		
		%if ( 201 == k )
		%	% missing this one!
		%	warning( "Missing 201!" );
		%	continue;
		%end
		
		vecX = load(sprintf( "%s/x%06d.txt", basePath, k ));
		vecF = load(sprintf( "%s/f%06d.txt", basePath, k ));
		assert( isreal(vecX) );
		assert( isreal(vecF) );
		assert( size(vecX,1) == sizeX );
		assert( size(vecF,1) == sizeF );
		assert( size(vecX,2) == 1 );
		assert( size(vecF,2) == 1 );
		%
		epsX = norm(vecX-vecX0);
		assert( epsX > 0.0 );
		rvecEpsX(k) = epsX;
		matJ(:,k) = (vecF-vecF0)/epsX;
		matY(:,k) = (vecX-vecX0)/epsX;
		vecFSumSq += vecF.^2;
	end
end
clear vecX;
clear vecF;
clear epsX;
clear k
clear m
clear n

%vecFSumSqMinScale = eps*max(vecFSumSq);
%vecFSumSq += vecFSumSqMinScale;
%vecFScale = sqrt(vecFSumSq/(sizeF*sizeK));
