CC=gcc
CFLAGS= -O2 -lm -lgsl -lgslcblas #-ggdb3
DEPS=
OBJS=main.o
OUT=main
GNUPLOT=gnuplot_2.sh

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)
	
.c.o:
	$(CC) -c -o $@ $< $(CFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	
run:
	./main

mem:
	valgrind ./main

graph:
	gnuplot $(GNUPLOT)

erase:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	rm -f *dat
	rm -f *txt
	rm -f *png

help:
	@echo -e "\033[1;33m============= makefile help =============\033[0m"
	@echo -e "\033[1;32mmake\033[0m       - kompiluje i linkuje pliki"
	@echo -e "\033[1;32mmake run\033[0m   - uruchamia program (plik wykonywalny)"
	@echo -e "\033[1;32mmake mem\033[0m   - sprawdza wycieki pamięci za pomocą valgrinda"
	@echo -e "\033[1;32mmake clean\033[0m - usuwa pliki .o oraz plik wykonywalny"
	@echo -e "\033[1;32mmake erase\033[0m - jak clean, ponadto: usuwa pliki tekstowe, pliki z danymi oraz wykresy"
	@echo -e "\033[1;32mmake graph\033[0m - wykonuje skrypt gnuplota i tworzy wykresy"
	@echo -e "\033[1;32mmake help\033[0m  - wyświetla pomoc"
	@echo -e "\033[1;33m=========================================\033[0m"