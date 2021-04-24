# GCCFLAGS = -pthread -c -O3 -flto
# GCCFLAGS = -std=c++11 -c -O3 -flto
GCCFLAGS = -std=c++17 -c -O3 -flto

evolve: evolve.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o
	g++ -o evolve evolve.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o -lpthread
demorun: demorun.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o
	g++ -o demorun demorun.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o -lpthread
singlerun: singlerun.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o
	g++ -o singlerun singlerun.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o Collide.o -lpthread
random.o: modules/random.cpp modules/random.h modules/VectorMatrix.h
	g++ $(GCCFLAGS) modules/random.cpp
TSearch.o: modules/TSearch.cpp modules/TSearch.h
	g++ $(GCCFLAGS) modules/TSearch.cpp
Worm.o: modules/Worm.cpp modules/Worm.h
	g++ $(GCCFLAGS) modules/Worm.cpp
Collide.o: modules/Collide.cpp modules/Collide.h
	g++ $(GCCFLAGS) modules/Collide.cpp
WormBody.o: modules/WormBody.cpp modules/WormBody.h modules/Collide.h
	g++ $(GCCFLAGS) modules/WormBody.cpp
NervousSystem.o: modules/NervousSystem.cpp modules/NervousSystem.h modules/VectorMatrix.h modules/random.h
	g++ $(GCCFLAGS) modules/NervousSystem.cpp
StretchReceptor.o: modules/StretchReceptor.cpp modules/StretchReceptor.h
	g++ $(GCCFLAGS) modules/StretchReceptor.cpp
Muscles.o: modules/Muscles.cpp modules/Muscles.h modules/VectorMatrix.h modules/random.h
	g++ $(GCCFLAGS) modules/Muscles.cpp
evolve.o: evolve.cpp modules/Worm.h modules/WormBody.h modules/StretchReceptor.h modules/Muscles.h modules/TSearch.h modules/Collide.h
	g++ $(GCCFLAGS) evolve.cpp
demorun.o: modules/Worm.h modules/WormBody.h modules/StretchReceptor.h modules/Muscles.h modules/TSearch.h modules/Collide.h
	g++ $(GCCFLAGS) demorun.cpp
singlerun.o: modules/Worm.h modules/WormBody.h modules/StretchReceptor.h modules/Muscles.h modules/TSearch.h modules/Collide.h
	g++ $(GCCFLAGS) singlerun.cpp
clean:
	rm *.o *.exe
wipedat:
	rm *.dat