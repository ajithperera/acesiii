./configure --enable-gnu \
            F77=mpif77 \
            FC=mpif77 \
            CC=mpicc \
            CXX=mpicxx \
            GNULIBS="-lgfortran -Wl,-Bstatic -lacml  -Wl,-Bdynamic" \
            GNUFLAGS="-L. -L/opt/amd/acml/5.2.0/gfortran64_fma4/lib" \
            FCFLAGS="-D__fortran -D__fortran77 -DMPIF2C -DMPI2 -Ofast -mavx -mfma4 -march=bdver1 -finit-local-zero -funroll-loops -falign-commons" \
            CFLAGS="-DMPIF2C -DC_SUFFIX -DMPI2 -DCB_SUFFIX -Ofast -mavx -mfma4 -march=bdver1" \
            CXXFLAGS="-DMPIF2C -DC_SUFFIX -DMPI2 -DCB_SUFFIX -Ofast -mavx -mfma4 -march=bdver1"


#FCFLAGS="-D__fortran -D__fortran77 -DMPIF2C -DMPI2 -O2 -pg -mavx -mfma4 -march=bdver1 -ffast-math -funroll-loops -fstack-arrays -finit-local-zero -falign-commons" \

#gfortran -o $@ -Wall -march=native -static -limf -ffast-math -funroll-
#loops -O3 $? 

#gfortran -O2 -funroll-loops --param max-unroll-times=4 -ftree-vectorize -ffast-math -march=corei7
#ifort -assume protect_parens -prec-div -prec-sqrt

#-O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -flto

#-ffast-math -funroll-loops -O3 -finline-limit=600 -fwhole-program -fstack-arrays -flto"
