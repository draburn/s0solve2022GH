function maxVal = mymax( varargin )
	for argIndex = 1: nargin
		x = varargin(argIndex){1};
		%
		xMax = max(x);
		n = 0;
		while (max(size(xMax))>1)
			xMax = max(xMax);
			n++;
			if (n>100)
				echo__xMax = xMax
				error("Something went wrong taking max() repeatedly.");
			end
		end
		%
		if ( 1==argIndex )
			maxVal = xMax;
		else
			maxVal = max([ xMax, maxVal ]);
		end
	end
	%
return;
end
