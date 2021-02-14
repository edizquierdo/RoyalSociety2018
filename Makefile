evolve: evolve.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o
	g++ -pthread -o evolve evolve.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o
random.o: modules/random.cpp modules/random.h modules/VectorMatrix.h
	g++ -c -O3 -flto modules/random.cpp
TSearch.o: modules/TSearch.cpp modules/TSearch.h
	g++ -c -O3 -flto modules/TSearch.cpp
Worm.o: modules/Worm.cpp modules/Worm.h
	g++ -c -O3 -flto modules/Worm.cpp
WormBody.o: modules/WormBody.cpp modules/WormBody.h
	g++ -c -O3 -flto modules/WormBody.cpp
NervousSystem.o: modules/NervousSystem.cpp modules/NervousSystem.h modules/VectorMatrix.h modules/random.h
	g++ -c -O3 -flto modules/NervousSystem.cpp
StretchReceptor.o: modules/StretchReceptor.cpp modules/StretchReceptor.h
	g++ -c -O3 -flto modules/StretchReceptor.cpp
Muscles.o: modules/Muscles.cpp modules/Muscles.h modules/VectorMatrix.h modules/random.h
	g++ -c -O3 -flto modules/Muscles.cpp
evolve.o: evolve.cpp modules/Worm.h modules/WormBody.h modules/StretchReceptor.h modules/Muscles.h modules/TSearch.h
	g++ -c -O3 -flto evolve.cpp
clean:
	rm *.o evolve
