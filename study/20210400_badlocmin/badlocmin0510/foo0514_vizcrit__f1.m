function f1 = foo0514_vizcrit__f1( x, y, pts )
	f1 = 5*x;
	%f1 = 10.0 + 0.2*(x.^2+y.^2);
	%f1 = x;
	%f1 = 100 + 0.1*( x.^2 + y.^2 );
	for n=1:size(pts,1)
		f1 += pts(n,4)*exp(-( (x-pts(n,1)).^2 + (y-pts(n,2)).^2 )/( 2.0*(pts(n,3).^2) ));
	end
return;
end
