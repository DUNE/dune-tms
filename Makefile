# The simplest makefile you've ever seen...
#
.PHONY: check-submodule all clean

default: all

# Check if git submodules need updating
check-submodule:
	@if git submodule status | egrep -q '^[-]|^[+]' ; then \
		echo "Will reinitialize git submodules"; \
	  	git submodule update --init; \
	fi

clean:
	make clean -C src
	make clean -C app

all: check-submodule
	make -j8 -C src
	make -j8 -C app
