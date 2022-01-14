function f_x_next = pickNextPoint_1D_0423( ...
  input__fa_x = [], ...
  input__ia_exceptionFlag = [], ...
  input__fa_f = [] )
  	% Should-be-precompiled...
	thisFile = "pickNextPoint_1D_0423.m";
	echo__input__fa_x = input__fa_x
	echo__input__fa_f = input__fa_f
	%
	% Local constants...
	myEps = abs(eps);
	mySqrtEps = sqrt(myEps);
	mySqrtSqrtEps = sqrt(mySqrtEps);
	assert( 0.0 < myEps );
	assert( 0.0 < mySqrtEps );
	assert( 0.0 < mySqrtSqrtEps );
	%
	% Indepdent input checks...
	assert(isrealarray(input__fa_x));
	assert(isrealarray(input__ia_exceptionFlag));
	assert(isrealarray(input__fa_f(0==input__ia_exceptionFlag)));
	%
	%size(input__fa_x)
	numPts = size(input__fa_x,2);
	% If no points, pick 0.0.
	if ( 0 == numPts )
		f_x_next = 0.0;
		return;
	end
	%
	% Cross-dependent inut checks...
	assert(isscalar(numPts));
	assert(isrealarray(input__fa_x,[1,numPts]));
	assert(isrealarray(input__ia_exceptionFlag,[1,numPts]));
	assert(isrealarray(input__fa_f(0==input__ia_exceptionFlag)));
	%
	% DO REAL WORK.
	% So, we have tried at least one point...
	assert(1<=numPts);
	%
	if ( 1 == numPts )
		% If only one point, require it to be valid,
		% and suggest another point.
		assert( 0 == input__ia_exceptionFlag(1) );
		if ( 0.0 < input__fa_x(1) )
			f_x_prev = input__fa_x(1);
			f_x_next = f_x_prev - abs(myEps*f_x_prev);
		elseif ( 0.0 > input__fa_x(1) )
			f_x_prev = input__fa_x(1);
			f_x_next = f_x_prev + abs(myEps*f_x_prev);
		else
			f_x_next = mySqrtSqrtEps;
		end
		return;
	end
	%
	%
	% Sort by x...
	[ sorted__fa_x, sorting__ia_n ] = sort( input__fa_x );
	sorted__ia_exceptionFlag = input__ia_exceptionFlag( sorting__ia_n );
	sorted__fa_f = input__fa_f( sorting__ia_n );
	echo__sorting__ia_n = sorting__ia_n
	echo__sorted__fa_x = sorted__fa_x
	echo__sorted__ia_exceptionFlag = sorted__ia_exceptionFlag
	echo__sorted__fa_f = sorted__fa_f
	%
	assert( 0 < sum(0==sorted__ia_exceptionFlag) );
	%
	trimIndexLo = 1;
	while ( 0 ~= sorted__ia_exceptionFlag(trimIndexLo) )
		trimIndexLo++;
	end
	if ( 1==trimIndexLo )
		haveLoBarrier = 0;
		xLoBarrier = -Inf;
	else
		haveLoBarrier = 1;
		xLoBarrier = sorted__fa_x(trimIndexLo-1);
	end
	%
	haveHiBarrier = 0;
	trimIndexHi = numPts;
	while ( 0 ~= sorted__ia_exceptionFlag(trimIndexHi) )
		trimIndexHi--;
	end
	if ( numPts==trimIndexHi )
		haveHiBarrier = 0;
		xHiBarrier = Inf;
	else
		haveHiBarrier = 1;
		xHiBarrier = sorted__fa_x(trimIndexHi+1);
	end
	%
	haveRangeBounds = [ haveLoBarrier, haveHiBarrier ]
	rangeBoundsVals = [ xLoBarrier, xHiBarrier ]
	trim__fa_x = sorted__fa_x(trimIndexLo:trimIndexHi)
	trim__fa_f = sorted__fa_f(trimIndexLo:trimIndexHi)
	trim__ia_exceptionFlag = sorted__ia_exceptionFlag(trimIndexLo:trimIndexHi)
	assert( 0 == sum(0~=trim__ia_exceptionFlag) );
	%
	f_x_next = pickNextPoint_1D_0423_trim( ...
	  trim__fa_x, trim__fa_f, haveRangeBounds, rangeBoundsVals );
return;
end

%!test
%!	assert( 0.0 == pickNextPoint_1D_0423 )

%!test
%!	assert( abs(pickNextPoint_1D_0423(0,0,0)-sqrt(sqrt(eps))) < sqrt(eps) );

%!test
%!	funch_ef = @(x)([ (0>=x)+(1<=x); x./((0<=x).*(1>=x)) ]);
%!	fa_x = linspace(-2,3,50);
%!	foobar = funch_ef(fa_x);
%!	ia_exceptionFlag = foobar(1,:);
%!	fa_f = foobar(2,:);

%!test
%!	funch_ef = @(x)([ (0>=x)+(1<=x); x./((0<=x).*(1>=x)) ]);
%!	fa_x = (3.0*rand(1,50))-1.0
%!	foobar = funch_ef(fa_x)
%!	ia_exceptionFlag = foobar(1,:)
%!	fa_f = foobar(2,:)
%!	
%!	pickNextPoint_1D_0423( fa_x, ia_exceptionFlag, fa_f );
