# The simplest makefile you've ever seen...
all:
	make -j8 -C src
	make -j8 -C app

clean:
	make clean -C src
	make clean -C app
