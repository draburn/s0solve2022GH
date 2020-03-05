function axpand( c = 0.2 )
	ax = axis();
	d = max([ ax(2)-ax(1), ax(4)-ax(3) ]);
	axis([ ...
	  ax(1) - c*d, ...
	  ax(2) + c*d, ...
	  ax(3) - c*d, ...
	  ax(4) + c*d ]);
return;
end
