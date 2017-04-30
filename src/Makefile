SOURCE=../
NAME=unimovib.exe

UNIMOVIB=$(SOURCE)$(NAME)

F90     = gfortran
#FFLAGS  = -O2 -w -static
FFLAGS  = -O2
LDFLAGS =

all:;
	make unimovibsrc

unimovibsrc:;
	make DIR=.             obj    FF=$(F90)
	make $(UNIMOVIB)


$(UNIMOVIB): objects/*.o
	rm -f $(UNIMOVIB)
	$(F90) $(FFLAGS) objects/*.o -o $(UNIMOVIB)

	@echo
	@echo Compiling -- $(UNIMOVIB) -- is done !
	@echo

obj:;
	ext=".f90";\
	for d in  $$DIR/*.f90; \
	do \
	  filename=`basename $$d $$ext`; \
	  make FILE=$$filename DIR=$$DIR FF=$$FF objects/$$filename.o; \
	done 

objects/$(FILE).o: $(DIR)/$(FILE).f90 ;
	$(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90;\
	mv $(FILE).o objects

clean:
	rm -f objects/*.o
	@echo delete objects/*.o
