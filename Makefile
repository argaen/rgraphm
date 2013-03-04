CXX=g++
CFLAGS= -I.
LIBS=-lgsl -lgslcblas -lm

all: recommendmake

recommendmake: main_recommender.cpp Node.cpp Link.cpp Group.cpp

		#$(CXX) -flto -O3 -o main_recommender $^ $(CFLAGS) $(LIBS)
		$(CXX) -g -o main_recommender $^ $(CFLAGS) $(LIBS)
