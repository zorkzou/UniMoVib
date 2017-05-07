SOURCE=../
NAME=unimovib.exe

UNIMOVIB=$(SOURCE)$(NAME)

F90     = gfortran
#FFLAGS  = -O2 -w -static
FFLAGS  = -O3
LDFLAGS =

all:;
	make unimovibsrc

unimovibsrc:;
	make DIR=.             obj90  FF=$(F90)
	make DIR=.             obj77  FF=$(F90)
	make $(UNIMOVIB)


$(UNIMOVIB): objects90/*.o objects77/*.o
	rm -f $(UNIMOVIB)
	$(F90) $(FFLAGS) objects90/*.o objects77/*.o -o $(UNIMOVIB)

	@echo
	@echo Compiling -- $(UNIMOVIB) -- is done !
	@echo

obj90:;
	ext=".f90";\
	for d in  $$DIR/*.f90; \
	do \
	  filename=`basename $$d $$ext`; \
	  make FILE=$$filename DIR=$$DIR FF=$$FF objects90/$$filename.o; \
	done 

obj77:;
	ext=".f";\
	for d in  $$DIR/*.f; \
	do \
	  filename=`basename $$d $$ext`; \
	  make FILE=$$filename DIR=$$DIR FF=$$FF objects77/$$filename.o; \
	done 

objects90/$(FILE).o: $(DIR)/$(FILE).f90 ;
	$(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90;\
	mv $(FILE).o objects90

objects77/$(FILE).o: $(DIR)/$(FILE).f ;
	$(FF) $(FFLAGS) -c $(DIR)/$(FILE).f;\
	mv $(FILE).o objects77

clean:
	rm -f objects90/*.o
	rm -f objects77/*.o
