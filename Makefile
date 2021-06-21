cython:
	export CFLAGS='-I$(shell python3 -c "import numpy;print(numpy.get_include())")' && python3 setup.py build_ext --inplace
	if [ -d "build" ]; then rm -Rf build; fi
	rm games/meangame.cpp

pytest:
	python3 -m pytest -q ./games/test/
