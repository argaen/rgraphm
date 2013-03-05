CXX=g++
CFLAGS= -I.
LIBS=-lgsl -lgslcblas -lm

all: recommendmake

recommendmake: main_recommender.cpp Node.cpp Link.cpp Group.cpp

		$(CXX) -Wall -flto -O3 -o main_recommender $^ $(CFLAGS) $(LIBS) #Production
		# $(CXX) -g -o main_recommender $^ $(CFLAGS) $(LIBS)    #Debugging compilation

alt: main.cpp Node.cpp Link.cpp Group.cpp
		$(CXX) -Wall -flto -O3 -o main_recommender $^ $(CFLAGS) $(LIBS) #Production

