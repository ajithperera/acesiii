#   
#   Main Makefile for sia
#

export
#   gmake all: Builds the xchssi executable.
CHSSI_EXE=xchssi

FILES:=$(wildcard *)
TARGET_DIRS:=$(shell for dir in $(FILES); \
                     do test -f $$dir/Makefile && echo $$dir; \
                     done)
#64BIT=1 # inherited from top-level make

all binclean libclean ppclean clean distclean: % : ;
	@for dir in $(TARGET_DIRS) ; \
	 do $(MAKE) -C $$dir $@ || exit 1 ; \
	 done

relink: binclean all

rebuild: libclean all

archive:


