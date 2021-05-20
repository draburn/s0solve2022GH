thisFile = "groot1d__getXNew";

numPts = size(xVals_raw,2);
assert( fevalCount == numPts );
if ( 0 == numPts )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'given - x1'." );
	xNew = x1;
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end
if ( 1 == numPts )
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'given - x2'." );
	xNew = x2;
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end

assert( 2 <= fevalCount );
xMin = min(xVals_raw);
xMax = max(xVals_raw);
fMin = min(fVals_raw);
fMax = max(fVals_raw);

% Handle "bounded".
if ( fMin * fMax < 0.0 )
	if ( 0 == boundedIter )
		msg_main( verbLev, thisFile, __LINE__, "Bounded a zero." );
	end
	n = 1;
	while ( fVals_sorted(n)*fVals_sorted(n+1) > 0.0 )
		n++;
	end
	assert( fVals_sorted(n)*fVals_sorted(n+1) < 0.0 );
	%
	if (0)
	%
	%
	% For simplicity, just alternate between lin and bisection.
	% More advanced methods would include using a cubic
	%  (after checking for good agreement with points further out),
	%  and using an "reduce distance between bounding points"
	%  that is smarter than just bisection.
	% But, all of this is irrelevant to higher dim anyway.
	boundedIter++; m = mod(boundedIter,2);
	if ( 0==m )
		msg_copious( verbLev, thisFile, __LINE__, "xNew via 'bounded - bisection'." );
		xNew = ( xVals_sorted(n) + xVals_sorted(n+1) )/2.0;
	else
		msg_copious( verbLev, thisFile, __LINE__, "xNew via 'bounded - linear interpolation'." );
		xNew = xVals_sorted(n) - fVals_sorted(n) ...
		 * ( xVals_sorted(n)-xVals_sorted(n+1) ) ...
		 / ( fVals_sorted(n)-fVals_sorted(n+1) );
	end
	%
	%
	else
	%
	%
	% Consider two quadratic models, look at "merit" of each.
	% Negative merit means "do not use".
	meritL = -1.0;
	meritR = -1.0;
	if (1)
	if ( 2 <= n )
		% Consider left model.
		x_opp = xVals_sorted(n+1); % Bound R.
		x_adj = xVals_sorted(n);   % Bound L.
		x_nbo = xVals_sorted(n-1); % Non-bounding.
		f_opp = fVals_sorted(n+1); % Bound R.
		f_adj = fVals_sorted(n);   % Bound L.
		f_nbo = fVals_sorted(n-1); % Non-bounding.
		%
		dxnorm = ( x_nbo - x_adj ) / ( x_adj - x_opp );
		dfnorm = ( f_nbo - f_adj ) / ( f_adj - f_opp );
		lx = log(dxnorm); % Prefer dxnorm ~ 1.
		lf = log(dfnorm/dxnorm); % Perfer slope ratio ~ 1.
		meritL = exp( -0.01*( lx^2 + lf^2 ) );
		msg_copious( verbLev, thisFile, __LINE__, ...
		  sprintf("dxnorm = %f.", dxnorm) );
		msg_copious( verbLev, thisFile, __LINE__, ...
		  sprintf("dfnorm = %f.", dfnorm) );
		msg_copious( verbLev, thisFile, __LINE__, ...
		  sprintf("meritL = %g.", meritL) );
		%
		if ( meritL < eps )
			meritL = -1.0;
		else
			vecX = [ x_opp; x_adj; x_nbo ]; % Note: decreasing in x.
			vecF = [ f_opp; f_adj; f_nbo ];
			matX = [ ones(3,1), vecX, vecX.^2 ];
			vecC = matX \ vecF;
			assert( isrealarray(vecC,[3,1]) );
			c0 = vecC(1);
			c1 = vecC(2);
			c2 = vecC(3);
			if ( abs(4.0*c0*c2) < sqrt(eps)*(c1^2) )
				% Use linear model.
				assert( 0.0 != c1 );
				xLeft = -c0/c1;
				assert( isrealscalar(xLeft) );
				msg_copious( verbLev, thisFile, __LINE__, ...
				  sprintf("xLeft = %f.", xLeft) );
				assert( xVals_sorted(n) < xLeft );
				assert( xLeft < xVals_sorted(n+1) );
			else
				discr = c1^2 - 4.0*c0*c2;
				assert( 0.0 <= discr );
				assert( 0.0 != c2 );
				xLeft = ( -c1 - sign(c1)*sqrt(discr) ) / ( 2.0 * c2 );
				assert( isrealscalar(xLeft) );
				msg_copious( verbLev, thisFile, __LINE__, ...
				  sprintf("xLeft = %f.", xLeft) );
				assert( xVals_sorted(n) < xLeft );
				assert( xLeft < xVals_sorted(n+1) );
			end
		end % End if ( meritL < 0.01 ) else.
	end % End if ( 2 <= n ).
	end
	if (1)
	if ( n+2 <= numPts )
		% Consider right model.
		x_opp = xVals_sorted(n);   % Bound L.
		x_adj = xVals_sorted(n+1); % Bound R.
		x_nbo = xVals_sorted(n+2); % Non-bounding.
		f_opp = fVals_sorted(n);   % Bound L.
		f_adj = fVals_sorted(n+1); % Bound R.
		f_nbo = fVals_sorted(n+2); % Non-bounding.
		%
		dxnorm = ( x_nbo - x_adj ) / ( x_adj - x_opp );
		dfnorm = ( f_nbo - f_adj ) / ( f_adj - f_opp );
		lx = log(dxnorm); % Prefer dxnorm ~ 1.
		lf = log(dfnorm/dxnorm); % Perfer slope ratio ~ 1.
		meritR = exp( -0.01*( lx^2 + lf^2 ) );
		%
		if ( meritR < eps )
			meritR = -1.0;
		else
			vecX = [ x_opp; x_adj; x_nbo ]; % Note: decreasing in x.
			vecF = [ f_opp; f_adj; f_nbo ];
			matX = [ ones(3,1), vecX, vecX.^2 ];
			vecC = matX \ vecF;
			assert( isrealarray(vecC,[3,1]) );
			c0 = vecC(1);
			c1 = vecC(2);
			c2 = vecC(3);
			if ( abs(4.0*c0*c2) < sqrt(eps)*(c1^2) )
				% Use linear model.
				assert( 0.0 != c1 );
				xRight = -c0/c1;
				assert( isrealscalar(xRight) );
				assert( xVals_sorted(n) < xRight );
				assert( xRight < xVals_sorted(n+1) );
			else
				discr = c1^2 - 4.0*c0*c2;
				assert( 0.0 <= discr );
				assert( 0.0 != c2 );
				xRight = ( -c1 - sign(c1)*sqrt(discr) ) / ( 2.0 * c2 );
				assert( isrealscalar(xRight) );
				assert( xVals_sorted(n) < xRight );
				assert( xRight < xVals_sorted(n+1) );
			end
		end % End if ( meritR < 0.01 ) else.
	end % End if ( n+2 <= numPts ).
	end
	%
	if ( meritL < 0.0 )
		if ( meritR < 0.0 )
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'bounded - linear interpolation'." );
			xNew = xVals_sorted(n) - fVals_sorted(n) ...
			 * ( xVals_sorted(n)-xVals_sorted(n+1) ) ...
			 / ( fVals_sorted(n)-fVals_sorted(n+1) );
		else
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'bounded - right-side quadratic interpolation'." );
			xNew = xRight;
		end
	else
		if ( meritR < 0.0 )
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'bounded - left-side quadratic interpolation'." );
			xNew = xLeft;
		else
			msg_copious( verbLev, thisFile, __LINE__, ...
			  "xNew via 'bounded - blended quadratic interpolation'." );
			xNew = ( meritL*xLeft + meritR*xRight ) / ( meritL + meritR );
		end
	end
	assert( isrealscalar(xNew) );
	assert( xVals_sorted(n) < xNew );
	assert( xNew < xVals_sorted(n+1) );
	%
	%
	end
	%
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end

% Handle typical case.
xNew = findGoodCand( xVals_sorted, fVals_sorted );
if ( isrealscalar(xNew) )
	xCapHi = xVals_sorted(end) + (xVals_sorted(end)-xVals_sorted(1));
	xCapLo = xVals_sorted(1)   - (xVals_sorted(end)-xVals_sorted(1));
	xNew = cap( xNew, xCapLo, xCapHi );
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'typical method - ???'." );
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end

% Handle "cylce: expand left, expand right, bisect"?
% Handle "all values are the same".
%%% HOW SHOULD THIS COMPETE WITH LARGER-SCALE TYPICAL CASE???
m = mod(desperationIter,5); desperationIter++;
if (0==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand right'." );
	xNew = xMax + 0.3*(xMax-xMin);
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
elseif (1==m)
	msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - expand left'." );
	xNew = xMin - 0.3*(xMax-xMin);
	thisFile = "RETURNING FROM groot1d__getXNew";
	return;
end
%
msg_copious( verbLev, thisFile, __LINE__, "xNew via 'desperation - bisect largest gap'." );
[ largestGapSize, largestGapIndex ] = max(diff(xVals_sorted));
xNew = (xVals_sorted(largestGapIndex) + xVals_sorted(largestGapIndex+1))/2.0;
thisFile = "RETURNING FROM groot1d__getXNew";
return;
