%  Function...
%    borderzoom( c = 0.2 )
%  Overview...
%    Part of myutil module.
%    A simple function to zoom in or out the current figure;
%    the zoom adds (or removes) a uniform border to the graph.
%  Input...
%    c: The factor by which to zoom in or out.
function borderzoom( c = 0.4 )
	assert(isrealscalar(c));
	ax = axis();
	d = 0.5 * max([ ax(2)-ax(1), ax(4)-ax(3) ]);
	axis([ ...
	  ax(1) - c*d, ...
	  ax(2) + c*d, ...
	  ax(3) - c*d, ...
	  ax(4) + c*d ]);
return;
end

%!test
%!	x = (0:100)/100;
%!	y = sin(2*pi*x);
%!	figure(gcf+1);
%!	plot(x,y,'o-');
%!	grid on;
%!	borderzoom;
