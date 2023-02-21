#import inspect
import linenumdemoB
from danutil import msg
#main_frame = inspect.currentframe()
#def msg(*arguments, **keywords):
#	print(f'[{__file__}.{inspect.stack()[1].lineno:05d}]', *arguments, **keywords)

msg('Hello from line 8 of linenumdemoA!')

linenumdemoB.dummyfunc()

msg('Hello from line 12 of linenumdemoA!')


def dummierfunc():
	msg('Hello from linenumdemoA.dummierfunc(), line 16 of linenumdemoB!')

msg('Hello from line 18 of linenumdemoA!')
dummierfunc()
