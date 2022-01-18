page_screen_output(0)
pkg load optim;
addpath( "myutil" );
addpath( "common" );
addpath( "numopt" );
commondefs;
%
%my_editor = "/usr/bin/codeblocks";
my_editor = "/usr/bin/mousepad";
EDITOR( [my_editor ' %s'] ); % Moot if not USING_GUI; see options within GUI.
%
switch 4
case 0
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0117lostinparadise/\".");
	addpath( "study/20211111hotCurves/0117lostinparadise/" );
case 1
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0101trine/\".");
	addpath( "study/20211111hotCurves/0101trine/" );
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0101trine/numopt0104/\".");
	addpath( "study/20211111hotCurves/0101trine/numopt0104/" );
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0101trine/testfunc2021/\".");
	addpath( "study/20211111hotCurves/0101trine/testfunc2021/" );
case 2
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0103onward/\".");
	addpath( "study/20211111hotCurves/0103onward/" );
case 3
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0110review/\".");
	addpath( "study/20211111hotCurves/0110review/" );
case 4
otherwise
	msg(__FILE__,__LINE__,"Invalid case.");
	error("Invalid case.");
end
