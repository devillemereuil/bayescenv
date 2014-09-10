# BayeScan makefile

bayescenv: start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 
	g++ -fopenmp -O3 -static -lpthread -o bayescenv start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 

start.o: start.cpp errors.cpp anyoption.h global_defs.h
	g++ -fopenmp -O3 -c start.cpp errors.cpp 

beta.o: beta.cpp global_defs.h
	g++ -fopenmp -O3 -c beta.cpp 
      
dirichlet.o: dirichlet.cpp global_defs.h
	g++ -fopenmp -O3 -c dirichlet.cpp 

RJupdates.o: RJupdates.cpp global_defs.h
	g++ -fopenmp -O3 -c RJupdates.cpp 

MHupdates.o: MHupdates.cpp global_defs.h
	g++ -fopenmp -O3 -c MHupdates.cpp 

likelihood.o: likelihood.cpp global_defs.h
	g++ -fopenmp -O3 -c likelihood.cpp 

read_write.o: read_write.cpp errors.cpp global_defs.h
	g++ -fopenmp -O3 -c read_write.cpp errors.cpp 

anyoption.o: anyoption.cpp anyoption.h 
	g++ -fopenmp -O0 -c anyoption.cpp 

clean: 
	rm *.o bayescenv