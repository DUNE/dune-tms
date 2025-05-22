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
	$(MAKE) clean -C src
	$(MAKE) clean -C app

all: check-submodule
	$(MAKE) -C src
	$(MAKE) -C app

sanitize: check-submodule
	$(MAKE) -C src sanitize
	$(MAKE) -C app sanitize
