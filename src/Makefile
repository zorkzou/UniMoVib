#
# Written by W. Zou <qcband@gmail.com>. Fri Nov 20 2020.
#

#DIR        = ../win32/
DIR        = ../
NAME       = unimovib.exe
EXE        = $(DIR)$(NAME)

OBJ        = object

FCOM       = gfortran
#            For GCC 10, -fallow-argument-mismatch is needed.
COMPFLAGS  = -O3
#            Looking for .o & .mod files here!
COMPFLAGS += -I$(OBJ) -J$(OBJ)

SRC90      = interface.f90  main.f90  rw.f90  symmetry.f90  thermo.f90  util.f90  vibfreq.f90
SRC77      = BLAS.f  LAPACK.f
MODS       = 

F90OBJ     = $(addprefix $(OBJ)/, $(SRC90:.f90=.o))
F77OBJ     = $(addprefix $(OBJ)/, $(SRC77:.f=.o))
MODOBJ     = $(addprefix $(OBJ)/, $(MODS:.f90=.o))

all: $(MODOBJ) $(F90OBJ) $(F77OBJ)
	@echo
	$(FCOM) $(COMPFLAGS) $(OBJ)/*.o -o $(EXE)
	@echo
	@echo "  --- $(NAME) has been compiled successfully!"
	@echo

$(MODOBJ): $(OBJ)/%.o: %.f90
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f90
#	mv *.mod $(OBJ)

$(F90OBJ): $(OBJ)/%.o: %.f90
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f90

$(F77OBJ): $(OBJ)/%.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

.PHONY: clean
clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.MOD
