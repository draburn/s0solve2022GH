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
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0205LandOfNod/\".");
	addpath( "study/20211111hotCurves/0205LandOfNod/" );
case 2
	msg(__FILE__,__LINE__,"Adding path \"study/20220208NoiseCompensation/0208OneDLinPlusRand/\".");
	addpath( "study/20220208NoiseCompensation/0208OneDLinPlusRand/" );
case 3
	msg(__FILE__,__LINE__,"Adding path \"study/20211111hotCurves/0210ThreeDim/\".");
	addpath( "study/20211111hotCurves/0210ThreeDim/" );
case 4
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0222alytJ/\".");
	addpath("study/20220222findLocMin/0222alytJ/");
otherwise
	msg(__FILE__,__LINE__,"Invalid case.");
	error("Invalid case.");
end
