import ctypes,os,platform,sys
from numpy.ctypeslib import ndpointer
from numpy import array

pyvernum = str(sys.version_info[0])+str(sys.version_info[1])
libpath = os.path.abspath(__file__)
libpath = os.path.realpath(libpath)
libpath = os.path.dirname(libpath)
libpath = os.path.dirname(libpath)

# try to figure out the name of the library
if platform.system() == "Darwin":
    libpath = os.path.join(libpath,"lzc.cpython-"+pyvernum+"m-darwin.so")
elif platform.system() == "Windows":
    libpath = os.path.join(libpath,"lzc.cp"+pyvernum+"-win_amd64.pyd")
else:
    libpath = os.path.join(libpath,"lzc.cpython-"+pyvernum+"m-x86_64-linux-gnu.so")

# path to DLL and argument types
_lzc = ctypes.CDLL(libpath)
_lzc.lz_complexity.argtypes = (ndpointer(ctypes.c_int), ctypes.c_int)
_lzc.lz_complexity2.argtypes = (ndpointer(ctypes.c_int), ctypes.c_int, ctypes.c_int)

def lz_complexity(s):
    global _lzc
    # coerce input arguments to ctypes
    s = array(s, dtype=ctypes.c_int)
    N = ctypes.c_int(len(s))
    lzc = _lzc.lz_complexity(s,N)
    return lzc

    """
    if platform.system() == 'Windows':
       #Load DLL file
       #lzcpythonDLL = ctypes.CDLL("<path to LZCINC.DLL>LZCINC.DLL")
       lzcpythonDLL = ctypes.CDLL("LZCINC.dll")
       #Make sure arguments are passable
       lzc = lzcpythonDLL.lz_complexity(p_s,N)
       del lzcpythonDLL

    #Default to Mac OS
    #Darwin working!!!
    else:
        libpath = os.path.abspath(__file__)
        libpath = os.path.realpath(libpath)
        libpath = os.path.dirname(libpath)
        libpath = os.path.join(libpath, "liblz_complexity.so")

        _lzcpython = ctypes.CDLL(libpath)
        _lzcpython.lz_complexity.argtypes= (ndpointer(ctypes.c_int), ctypes.c_int)
        _lzcpython.lz_complexity.restype = ctypes.c_int

        #Library call
        lzc = _lzcpython.lz_complexity(s, N, res)

    return lzc
    """
def lz_complexity2(s, threshold):
    global _lzc
    # coerce input arguments to ctypes
    s = array(s, dtype=ctypes.c_int)
    N = ctypes.c_int(len(s))
    lzc = _lzc.lz_complexity2(s,N,ctypes.c_int(threshold))
    return lzc