cython:
	export CFLAGS='-I$(shell python -c "import numpy;print(numpy.get_include())")' && python setup.py build_ext --inplace
	if [ -d "build" ]; then rm -Rf build; fi
	rm games/meangame.cpp

pytest: spatgames.*.so
	python -m pytest -q ./games/test/
