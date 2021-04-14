from distutils.core import Extension, setup

extension_mod = Extension("mod_p9m3upv6",
		[ r'mod_p9m3upv6_wrapper.c' ],
		extra_objects = [r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__/bind_c_mod_p9m3upv6.o',
				r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__/mod_p9m3upv6.o'],
		include_dirs = [r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__',
				r'/home/cherif/anaconda3/lib/python3.8/site-packages/numpy/core/include'],
		libraries = [r'gfortran'],
		library_dirs = [r'/usr/lib/gcc/x86_64-linux-gnu/9'],
		extra_compile_args = [r'-O3',
				r'-fPIC',
				r'-I"/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__"'],
		extra_link_args = [])

setup(name = "mod_p9m3upv6", ext_modules=[extension_mod])