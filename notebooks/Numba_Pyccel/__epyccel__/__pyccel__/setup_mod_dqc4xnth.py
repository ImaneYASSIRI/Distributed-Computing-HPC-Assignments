from distutils.core import Extension, setup

extension_mod = Extension("mod_dqc4xnth",
		[ r'mod_dqc4xnth_wrapper.c' ],
		extra_objects = [r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__/bind_c_mod_dqc4xnth.o',
				r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__/mod_dqc4xnth.o'],
		include_dirs = [r'/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__',
				r'/home/cherif/anaconda3/lib/python3.8/site-packages/numpy/core/include'],
		libraries = [r'gfortran'],
		library_dirs = [r'/usr/lib/gcc/x86_64-linux-gnu/9'],
		extra_compile_args = [r'-O3',
				r'-fPIC',
				r'-I"/home/cherif/Bureau/HPC/Distributed-Computing-HPC-Assignments/notebooks/Numba_Pyccel/__epyccel__/__pyccel__"'],
		extra_link_args = [])

setup(name = "mod_dqc4xnth", ext_modules=[extension_mod])