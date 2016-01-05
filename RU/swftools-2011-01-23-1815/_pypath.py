import distutils
import distutils.sysconfig
import sys
sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=0,standard_lib=0))
