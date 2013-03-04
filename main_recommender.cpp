#include <iostream>
#include <string.h>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <boost/unordered_map.hpp>
#include <vector>

#include "Node.h"
#include "Link.h"
#include "Group.h"

#define STEPS 100

typedef boost::unordered_map<int, Node> Hash_Map;
//typedef std::vector<Node> Group;
typedef boost::unordered_map<int, Group> Groups;

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

/* ##################################### */
double Group2GroupH(Group *g1, Group *g2, int K){

	double H = 0.0;

	H += gsl_sf_lnfact(g1->g2glinks[g2->get_id()][0] + K -1);
	for(int i=1; i<K+1; ++i){
		H -= gsl_sf_lnfact(g1->g2glinks[g2->get_id()][i]);
	}
    return H;
}


/* ##################################### */
double HKState(int K, Groups *g1, Groups *g2, Hash_Map d1, Hash_Map d2){
	
	
	double H = 0.0;
	int nnod1=0, nnod2=0;

	for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
		for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
			H += Group2GroupH(&(it1->second), &(it2->second), K);
		
	H -= gsl_sf_lnfact(d1.size() - g1->size()) + gsl_sf_lnfact(d2.size() - g2->size());
	
	return H;

}


/* ##################################### */
void CreateRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, gsl_rng *rgen, int K, bool random){

	int nweight[K+1];
	memset(nweight, 0, sizeof(int)*(K+1));
	int nlinks = 0, ngrouplinks = 0;
	int group;
	
	for (int i = 1; i<(*d1).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d1).size())) + 1) : group=i;
		(*d1)[i].set_group(group);
		((*groups1)[group]) = Group(group, K);
		(*groups1)[group].members[group]=(*d1)[i];
	}
	
	for (int i = 1; i<(*d2).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d2).size())) + 1) : group=i;
		(*d2)[i].set_group(group);
		((*groups2)[group]) = Group(group, K);
		(*groups2)[group].members[group]=(*d2)[i];
	}

	for (Groups::iterator g1 = groups1->begin(); g1 != groups1->end(); ++g1){
        for (Groups::iterator g2 = groups2->begin(); g2 != groups2->end(); ++g2){

		//Groups 2 info filling
			for(GroupNodes::iterator it1 = g2->second.members.begin(); it1 != g2->second.members.end(); ++it1)
    		    for (Links::iterator nit = it1->second.neighbours.begin(); nit != it1->second.neighbours.end(); ++nit){
        		    ++nlinks;
            		for(GroupNodes::iterator it2 = g1->second.members.begin(); it2 != g1->second.members.end(); ++it2)
                		if (nit->second.get_id() == it2->second.get_id())
                    		++ngrouplinks;
                
	        	}

			g2->second.set_nlinks(nlinks);
		    g2->second.g2glinks[g1->second.get_id()][0] = ngrouplinks;
		    for(int i=1; i<K+1; ++i){
		        g2->second.g2glinks[g1->second.get_id()][i] = nweight[i];
		        g2->second.kratings[i] += nweight[i];
		    }

		    nlinks = ngrouplinks = 0;
		    memset( nweight, 0, (K+1)*sizeof(int) );

		//Groups 1 info filling

			for(GroupNodes::iterator it1 = g1->second.members.begin(); it1 != g1->second.members.end(); ++it1)
		        for (Links::iterator nit = it1->second.neighbours.begin(); nit != it1->second.neighbours.end(); ++nit){
        		    ++nlinks;
		            for(GroupNodes::iterator it2 = g2->second.members.begin(); it2 != g2->second.members.end(); ++it2)
		                if (nit->second.get_id() == it2->second.get_id())
        		            ++ngrouplinks;
                
        		}

		    g1->second.set_nlinks(nlinks);
		    g1->second.g2glinks[g2->second.get_id()][0] = ngrouplinks;
		    for(int i=1; i<K+1; ++i){
		        g1->second.g2glinks[g2->second.get_id()][i] = nweight[i];
		        g1->second.kratings[i] += nweight[i];
		    }
		}
	}
}

/* ##################################### */
int MCStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, double *H, int K){

	bool visitedgroup[K+1];
	memset( visitedgroup, false, (K+1)*sizeof(bool) );
	Group *src_g, *dest_g;
	int newgrp, oldgrp;
	int dice;
	int id;
	Node n;
	double dH = 0.0;
	int nr;

	dice = floor(gsl_rng_uniform(stepgen) * d1->size() + 1);
	n = d1->at(dice);
	oldgrp = n.get_group();
	newgrp = floor(gsl_rng_uniform(groupgen) * d1->size() + 1);
	while ( newgrp == oldgrp ) 
		newgrp = floor(gsl_rng_uniform(groupgen) * d1->size() + 1);

	src_g = &(*g1)[oldgrp];
	dest_g = &(*g1)[newgrp];

	for (Links::iterator it = n.neighbours.begin(); it != n.neighbours.end(); ++it){
		if (!visitedgroup[it->second.get_id()]){
		    id = (*d2)[it->second.get_id()].get_group();
			dH -= gsl_sf_lnfact(src_g->g2glinks[id][0] + K - 1);
			dH -= gsl_sf_lnfact(dest_g->g2glinks[id][0] + K -1);
			for (int i = 1; i<K+1; ++i){
				nr = src_g->g2glinks[id][i];
				dH -= gsl_sf_lnfact(nr);
				nr = dest_g->g2glinks[id][i];
				dH -= gsl_sf_lnfact(nr);		
			}
			visitedgroup[it->second.get_id()] = true;
		}
	}

	//if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
	//	dH += gsl_sf_lnfact(d1->size() - g1->size());

	//src_g->remove_node(n);
	//dest_g->add_node(n);
	(*d1)[n.get_id()].set_group(dest_g->get_id());

	for (Links::iterator it = n.neighbours.begin(); it != n.neighbours.end(); ++it){
		if (visitedgroup[it->second.get_id()]){
	        id = (*d2)[it->second.get_id()].get_group();
    	    dH += gsl_sf_lnfact(src_g->g2glinks[id][0] + K - 1);
	        dH += gsl_sf_lnfact(dest_g->g2glinks[id][0] + K -1);
    	    for (int i = 1; i<K+1; ++i){
	            nr = src_g->g2glinks[id][i];
    	        dH += gsl_sf_lnfact(nr);
        	    nr = dest_g->g2glinks[id][i];
            	dH += gsl_sf_lnfact(nr);
        	}
			visitedgroup[it->second.get_id()] = false;
		}
    }	

	//if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
       // dH -= gsl_sf_lnfact(d1->size() - g1->size());

	//Test if the movement is useful for the system
	if ( dH <= 0.0 || gsl_rng_uniform(stepgen) < exp(-dH)){
//		std::cout << "Next dH: "<< dH <<"\n";	
		*H += dH;
	}
	//else undo the movement
	else{
//		std::cout << "Movement not accepted: " << dH << "\n";
//		dest_g->remove_node(n);
//        src_g->add_node(n);
   		(*d1)[n.get_id()].set_group(src_g->get_id());
	}

	return 0.0;
}


/* ##################################### */
int main(int argc, char **argv){

	int id1, id2, weight;
	double H;
	Hash_Map d1, d2;
	Hash_Map d1c, d2c;
	Groups groups2, groups1;
	gsl_rng *rand_groups, *rand_step;
	char* tFileName;
	char* qFileName;
	int groupseed, stepseed;
	int mark;

	parseArguments(argc, argv, &tFileName, &qFileName, &stepseed, &groupseed, &mark);
	rand_step = gsl_rng_alloc(gsl_rng_mt19937);
	rand_groups = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rand_step, stepseed);
	gsl_rng_set(rand_groups, groupseed);

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


	d1c = d1;
	d2c = d2;	

	CreateRandomGroups(&d1c, &d2c, &groups1, &groups2, rand_groups, mark, false);
	H = HKState(mark, &groups1, &groups2, d1c, d2c);
	std::cout << "Initial H: "<< H <<"\n";
	for(int i=0; i<STEPS; i++){
//		MCStepKState(&groups1, &groups2, &d1c, &d2c, rand_step, rand_groups, &H, mark);
		std::cout << H << "    " << HKState(mark, &groups1, &groups2, d1c, d2c) << "\n";
	}

	std::cout << " ################ G R O U P S 1 ############### \n";
    for (Groups::iterator it = (groups1).begin(); it != (groups1).end(); ++it){
		std::cout << "[Group:" << it->second.get_id();
        for(GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << ", Node:" << it1->second.get_id() << " Links: ";
                for (Links::iterator nit = it1->second.neighbours.begin(); nit != it1->second.neighbours.end(); ++nit)
                    std::cout << nit->second.get_id() << ", ";

        }
		std::cout << "] \nTotalLinks: " <<it->second.get_nlinks() << "\n";
		std::cout << "Kratings: ";
		for (int i=1; i<mark+1; ++i)
			std::cout << it->second.kratings[i] << " ";
		std::cout << "\n";
			
	}

    /*std::cout << " ################ G R O U P S 2 ################ \n";
    for (Groups::iterator it = (groups2).begin(); it != (groups2).end(); ++it){

        for(GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << "[ Group:" << it->second.get_id() << ", Node:" << (*it1).get_id() << " Links: ";
                 for (Links::iterator nit = (*it1).neighbours.begin(); nit != (*it1).neighbours.end(); ++nit)
                    std::cout << nit->second.get_id() << ", ";
                std::cout << "] \n";
        }
		std::cout << "TotalLinks: " <<it->second.get_nlinks() << "\n";
        std::cout << "Kratings: ";
        for (int i=1; i<mark+1; ++i)
            std::cout << it->second.kratings[i] << " ";
        std::cout << "\n";
	}*/

    std::cout << " ################################################# \n";

	std::cout << "Next H: "<< H <<"\n";

	gsl_rng_free(rand_step);
	gsl_rng_free(rand_groups);
	return 0;

}