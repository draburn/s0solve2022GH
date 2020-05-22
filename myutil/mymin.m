function minVal = mymin( varargin )
	for argIndex = 1: nargin
		x = varargin(argIndex){1};
		xMin = min(x);
		n = 0;
		while (max(size(xMin))>1)
			xMin = min(xMin);
			n++;
			if (n>100)
				echo__xMin = xMin
				error("Something went wrong taking min() repeatedly.");
			end
		end
		%
		if ( 1==argIndex )
			minVal = xMin;
		else
			minVal = min([ xMin, minVal ]);
		end
	end
	%
return;
end
