CXX=g++
CFLAGS= -Wall -O3 -std=c++11 -I.
# CFLAGS= -Wall -g -O3 -flto -I.
LIBS=-lgsl -lgslcblas -lm

all: production

debug: recommender.cpp Node.cpp Link.cpp Group.cpp utils.cpp

		# $(CXX) $(CFLAGS) -o recommender $^ $(LIBS) #Production
		$(CXX) -g -o recommender $^ $(CFLAGS) $(LIBS)    #Debugging compilation

production: recommender.cpp Node.cpp Link.cpp Group.cpp utils.cpp
		$(CXX) -Wall -flto -O3 -o recommender $^ $(CFLAGS) $(LIBS) #Production

