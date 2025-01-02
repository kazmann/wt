CC=/usr/bin/gcc
I=../include
LIB=../lib

wt2:	new_tephra.o parameters_strat_mac.h tephra_calc.o windy.o
		gcc -o wt2 new_tephra.o tephra_calc.o windy.o

new_tephra.o:	new_tephra.c common_structures_strat.h
		gcc -c new_tephra.c

tephra_calc.o:		tephra_calc.c common_structures_strat.h
		gcc -c tephra_calc.c
		
windy.o:	windy.c parameters_strat_mac.h common_structures_strat.h prototypes_strat.h
		gcc -c windy.c

clean:
		rm -fv *.o *.bak wt2