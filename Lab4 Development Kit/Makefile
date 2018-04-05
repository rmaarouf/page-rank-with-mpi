all: datatrim serialtester main

datatrim: datatrim.c
	gcc datatrim.c -o datatrim

serialtester: serialtester.c
	gcc serialtester.c Lab4_IO.c -o serialtester -lm

main: main.c
	mpicc main.c Lab4_IO.c -o main -lm

clean:
	rm -f datatrim
	rm -f main
	rm -f data_input
	rm -f data_output
	rm -f serialtester
