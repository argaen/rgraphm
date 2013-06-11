#include <iostream>
#include <string.h>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <iomanip>
#include <sys/time.h>


#include "Node.h"
#include "Link.h"
#include "Group.h"

#define STEPS 1000

typedef boost::unordered_map<int, double> LnFactList;
typedef boost::unordered_map<int, Node> Hash_Map;
//typedef std::vector<Node> Group;
typedef boost::unordered_map<int, Group> Groups;
typedef boost::unordered_map<int, Group> Groups;
typedef boost::multi_array<int, 3> GGLinks;
typedef GGLinks::index GG_index;

/* ##################################### */
void parseArguments(int argc, char **argv, char** inFile, char** qFile, int* stepseed, int* groupseed, int* mark) {
	if (argc != 11) {

		fprintf (stderr, "Usage: main_recommender -q queryFile -t trainFile -s stepseed -g groupsseed -m mark \n");
        exit(1);
	}

	int c;
	while ((c = getopt (argc, argv, "t:q:s:g:m:")) != -1)
    	switch (c) {
			case 't':
           		*inFile = optarg;
	            break;
			case 'q':
    	        *qFile = optarg;
        	    break;
			case 'g':
            	*groupseed = atoi(optarg);
             	break;
			case 's':
				*stepseed = atoi(optarg);
				break;
			case 'm':
				*mark = atoi(optarg);
				break;
			case '?':
				if (optopt == 't' || optopt == 'q' || optopt == 'm' || optopt == 's' || optopt == 'g')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);

				exit(1);

			default:

				fprintf (stderr, "Usage: main_recommender -q queryFile -t trainFile -s seed -m mark \n");
				exit(1);
		}
}

int GetDecorrelationKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups){

    int nrep = 10, step;
    int x1, x2;

    x2 = (nnod1 + nnod2) / 5;
    if (x2 < 10)
        x2 = 10;
    x1 = x2 / 4;

    for (int i=0; i<nrep; i++){

        printf("Estimating decorrelation time %d/%d", i,nrep);

        for (step=0; step<=x2; step++){
            MCStepKState(groups1, groups2, d1c, d2c, step_randomizer, groups_randomizer, H, mark, lnfactlist, logsize, nnod1, nnod2, gglinks);

        }


    }


}


/* ##################################### */
double Group2GroupH(Group *g1, Group *g2, int K, GGLinks *gglinks){
	double H = 0.0;

	/* H += gsl_sf_lnfact(g1->g2glinks[g2->get_id()][0] + K -1); */
	H += gsl_sf_lnfact((*gglinks)[g1->get_id()][g2->get_id()][0] + K -1);
	for(int i=1; i<K+1; ++i){
		/* H -= gsl_sf_lnfact(g1->g2glinks[g2->get_id()][i]); */
        H -= gsl_sf_lnfact((*gglinks)[g1->get_id()][g2->get_id()][i]);
	}
    return H;
}


/* ##################################### */
double HKState(int K, Groups *g1, Groups *g2, int d1_size, int d2_size, GGLinks *gglinks){
	double H = 0.0;

	for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
		for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
			H += Group2GroupH(&(it1->second), &(it2->second), K, gglinks);
		
	H -= gsl_sf_lnfact(d1_size - g1->size()) + gsl_sf_lnfact(d2_size - g2->size());
	
	return H;

}


/* ##################################### */
void CreateRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, gsl_rng *rgen, int K, bool random, GGLinks *gglinks){
	int nweight[K+1];
	int ngrouplinks;
	int group;

	for (int i = 1; i<(*d1).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d1).size())) + 1) : group=i;
		(*d1)[i].set_group(group);
		((*groups1)[group]) = Group(group, K);
		(*groups1)[group].members[i]=&((*d1)[i]);
	}
	
	for (int i = 1; i<(*d2).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d2).size())) + 1) : group=i;
		(*d2)[i].set_group(group);
		((*groups2)[group]) = Group(group, K);
		(*groups2)[group].members[i]=&((*d2)[i]);
	}

	for (Groups::iterator g1 = groups1->begin(); g1 != groups1->end(); ++g1){
        for (Groups::iterator g2 = groups2->begin(); g2 != groups2->end(); ++g2){
        memset(nweight, 0, sizeof(int)*(K+1));
        ngrouplinks = 0;

		//Groups 1 info filling

			for(GroupNodes::iterator it1 = g1->second.members.begin(); it1 != g1->second.members.end(); ++it1)
		        for (Links::iterator nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); ++nit){
		            for(GroupNodes::iterator it2 = g2->second.members.begin(); it2 != g2->second.members.end(); ++it2)
		                if (nit->second.get_id() == it2->second->get_id()){
        		            ++ngrouplinks;
                            nweight[nit->second.get_weight()]++;
                        }
        		}

            (*gglinks)[g1->second.get_id()][g2->second.get_id()][0] = ngrouplinks;
		    for(int i=1; i<K+1; ++i)
                (*gglinks)[g1->second.get_id()][g2->second.get_id()][i] = nweight[i];
		}
	}
}

/* ##################################### */

double* GenLogFactList(int size){
    double* LogFactList = (double*) calloc(size, sizeof(double));

    for (int i = 0; i<size; i++)
        LogFactList[i] = gsl_sf_lnfact(i);

    return LogFactList;
}

double LogFact(int key, int size, double* LogFactList){
    if (size<key)
        return gsl_sf_lnfact(key);
    else
        return LogFactList[key];
}

int MCStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinks){


    struct timeval stop, start;
    gettimeofday(&start, NULL);
    Groups *g;
    Hash_Map *d_move, *d_nomove;
    int bnnod = (nnod1>nnod2) ? nnod1+1 : nnod2+1;

	bool visitedgroup[bnnod]; 
    for (int i=0; i < bnnod; i++)
        visitedgroup[i]=false;
	Group *src_g, *dest_g;
	int newgrp, oldgrp, dice, set_size_move, id;
    bool set_ind;
	Node *n;
	double dH = 0.0;


    if (gsl_rng_uniform(stepgen) < (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2)){
        g = g1;d_move = d1;d_nomove = d2;set_ind = true;
    }else{
        g = g2;d_move = d2;d_nomove = d1;set_ind = false;
    }
    set_size_move = d_move->size();


    dice = gsl_rng_uniform(stepgen) * set_size_move + 1;
	n = &(d_move->at(dice));
	oldgrp = n->get_group();
	newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;
	while ( newgrp == oldgrp ) 
		newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;

	src_g = &(*g)[oldgrp];
	dest_g = &(*g)[newgrp];

	for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        id = (*d_nomove)[it->second.get_id()].get_group();
		if (!visitedgroup[id]){
            if (set_ind){
                dH -= LogFact((*gglinks)[src_g->get_id()][id][0] + K - 1, logsize, lnfactlist);
                dH -= LogFact((*gglinks)[dest_g->get_id()][id][0] + K - 1, logsize, lnfactlist);
            }else{
                dH -= LogFact((*gglinks)[id][src_g->get_id()][0] + K - 1, logsize, lnfactlist);
                dH -= LogFact((*gglinks)[id][dest_g->get_id()][0] + K - 1, logsize, lnfactlist);
            }

			for (int i = 1; i<K+1; ++i){
                if (set_ind){
                    dH -= -LogFact((*gglinks)[src_g->get_id()][id][i], logsize, lnfactlist);
                    dH -= -LogFact((*gglinks)[dest_g->get_id()][id][i], logsize, lnfactlist);
                }else{
                    dH -= -LogFact((*gglinks)[id][src_g->get_id()][i], logsize, lnfactlist);
                    dH -= -LogFact((*gglinks)[id][dest_g->get_id()][i], logsize, lnfactlist);
                }
			}
			visitedgroup[id] = true;
		}
	}

	if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
		dH += LogFact(set_size_move - g1->size(), logsize, lnfactlist);

    if (set_ind){
        src_g->remove_node_s1(n, d_nomove, gglinks);
        dest_g->add_node_s1(n, d_nomove, gglinks);
    }else{
        src_g->remove_node_s2(n, d_nomove, gglinks);
        dest_g->add_node_s2(n, d_nomove, gglinks);
    }

	for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        id = (*d_nomove)[it->second.get_id()].get_group();
		if (visitedgroup[id]){
            if (set_ind){
                dH += LogFact((*gglinks)[src_g->get_id()][id][0] + K - 1, logsize, lnfactlist);
                dH += LogFact((*gglinks)[dest_g->get_id()][id][0] + K - 1, logsize, lnfactlist);
            }else{
                dH += LogFact((*gglinks)[id][src_g->get_id()][0] + K - 1, logsize, lnfactlist);
                dH += LogFact((*gglinks)[id][dest_g->get_id()][0] + K - 1, logsize, lnfactlist);
            }

			for (int i = 1; i<K+1; ++i){
                if (set_ind){
                    dH += -LogFact((*gglinks)[src_g->get_id()][id][i], logsize, lnfactlist);
                    dH += -LogFact((*gglinks)[dest_g->get_id()][id][i], logsize, lnfactlist);
                }else{
                    dH += -LogFact((*gglinks)[id][src_g->get_id()][i], logsize, lnfactlist);
                    dH += -LogFact((*gglinks)[id][dest_g->get_id()][i], logsize, lnfactlist);
                }
			}
			visitedgroup[id] = false;
		}
    }	

	if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
		dH -= LogFact(set_size_move - g1->size(), logsize, lnfactlist);

	if ( dH <= 0.0 || gsl_rng_uniform(stepgen) < exp(-dH)){
		*H += dH;
	}else{
        if(set_ind){
            dest_g->remove_node_s1(n, d_nomove, gglinks);
            src_g->add_node_s1(n, d_nomove, gglinks);
        }else{
            dest_g->remove_node_s2(n, d_nomove, gglinks);
            src_g->add_node_s2(n, d_nomove, gglinks);
        }
	}
    gettimeofday(&stop, NULL);
    printf("Time %lu\n", stop.tv_usec - start.tv_usec);


	return 0;
}

void printGroups(Groups g, int mark){
    for (Groups::iterator it = (g).begin(); it != (g).end(); ++it){
        std::cout << "[Group:" << it->second.get_id();
        for(GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << ", Node:" << it1->second->get_id() << " Links: ";
                for (Links::iterator nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); ++nit)
                    std::cout << nit->second.get_id() << ", ";

        }
        std::cout << "\n";
    }
}


/* ##################################### */
int main(int argc, char **argv){

    int decorStep;
	int id1, id2, weight, logsize=100;
	double H;
    double *lnfactlist;
	Hash_Map d1, d2;
	Hash_Map d1c, d2c;
	Groups groups2, groups1;
	gsl_rng *groups_randomizer, *step_randomizer;
	char* tFileName;
	char* qFileName;
	int groupseed, stepseed;
	int mark;
    int nnod1, nnod2;

	parseArguments(argc, argv, &tFileName, &qFileName, &stepseed, &groupseed, &mark);

    std::cout << stepseed << " " << tFileName << " " << qFileName << " " << groupseed << " " << mark << "\n";

    
	step_randomizer = gsl_rng_alloc(gsl_rng_mt19937);
	groups_randomizer = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(step_randomizer, stepseed);
	gsl_rng_set(groups_randomizer, groupseed);

	std::ifstream inFile(tFileName);
	if (inFile.is_open()){
		while (inFile >> id1 >> id2 >> weight){
			d1.insert(Hash_Map::value_type(id1,Node(id1)));
			d1[id1].neighbours.insert(Links::value_type(id2,Link(id2, weight)));
			d2.insert(Hash_Map::value_type(id2,Node(id2)));
			d2[id2].neighbours.insert(Links::value_type(id1,Link(id1, weight)));
			

		}
		inFile.close();

	}else
		std::cout << "Couldn't open file " << tFileName << "\n";
    
    lnfactlist = GenLogFactList(logsize);


	d1c = d1;
	d2c = d2;	
    nnod1 = d1.size();
    nnod2 = d2.size();

    GGLinks gglinks(boost::extents[nnod1+2][nnod2+2][mark+1]);


    GetDecorrelationKState();

	CreateRandomGroups(&d1c, &d2c, &groups1, &groups2, groups_randomizer, mark, false, &gglinks);


	H = HKState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks);
    double TH;
	std::cout << "Initial H: "<< H <<"\n";
	for(int i=0; i<STEPS; i++){
		MCStepKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks);
        TH = HKState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks);
        /* std::cout << std::setprecision(20) << H << "    " << "\n"; */
        std::cout << std::setprecision(20) << H << "    " << TH << "\n";
        /* printf("############GROUPS1################\n"); */
        /* printGroups(groups1, mark); */
        /* printf("############GROUPS2################\n"); */
        /* printGroups(groups2, mark); */
        /* printf("###################################\n"); */
	}

    free(lnfactlist);
	gsl_rng_free(step_randomizer);
	gsl_rng_free(groups_randomizer);
	return 0;

}

