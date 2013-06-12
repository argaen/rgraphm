#include <iostream>
#include <string.h>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multiroots.h>
#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>
#include <vector>
#include <iomanip>
#include <sys/time.h>


#include "Node.h"
#include "Link.h"
#include "Group.h"

#define STEPS 10

typedef boost::unordered_map<int, double> LnFactList;
typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map<int, Group> Groups;
typedef boost::multi_array<int, 3> GGLinks;

struct pair
{
  double x0;
  double y0;
  double x1;
  double y1;
};

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


int
ExponentialRootF(const gsl_vector *params, void *points, gsl_vector *f) {
    const double a = gsl_vector_get(params, 0);
    const double b = gsl_vector_get(params, 1);

    double x0 = ((struct pair *)points)->x0;
    double y0 = ((struct pair *)points)->y0;
    double x1 = ((struct pair *)points)->x1;
    double y1 = ((struct pair *)points)->y1;

    const double r0 = y0 - a - (1. - a) * exp(-x0 / b);
    const double r1 = y1 - a - (1. - a) * exp(-x1 / b);

    gsl_vector_set(f, 0, r0);
    gsl_vector_set(f, 1, r1);

    return GSL_SUCCESS;
}

double getDecay(int nnod, double x1, double x2, double y1, double y2) {

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int status;
    size_t i, iter = 0;
    const size_t n = 2;
    struct pair p = {x1, y1, x2, y2};
    gsl_multiroot_function f = {&ExponentialRootF, n, &p};
    double x_init[2] = {y2, sqrt(nnod)};
    gsl_vector *x = gsl_vector_alloc(n);
    double result;

    for (i=0; i<n; i++)
        gsl_vector_set(x, i, x_init[i]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);

    do {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      if (status)   /* check if solver is stuck */
        break;
      status =gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);

    if (strcmp(gsl_strerror(status), "success") != 0)
        result = -1;
    else
        result = gsl_vector_get(s->x, 1);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return result;
}

double mean(double* data, int N) {
    double acc = 0.0;
    for (int i=0; i<N; ++i)
        acc+=data[i];

    return acc/N;
}

double stddev(double *data, int N) {
    int i;
    double m = 0.0, m2 = 0.0;
    double s;

    for (i=0; i<N; i++) {
        m += data[i];
        m2 += data[i] * data[i];
    }

    m /= (double)(N);
    m2 /= (double)(N);
    if (m2 - m * m > 0.0)
        s = sqrt(m2 - m * m);
    else
        s = 0.0;

    return s;
}

double mutualInfo(Groups *g, Groups *gt, int nnod1, int nnod2){
    int S = 0, S1 = 0, S2 = 0, S12 = 0;
    double H1 = 0.0, H2 = 0.0, H12 = 0.0;

    S = nnod1;

    for (Groups::iterator it = (*g).begin(); it != (*g).end(); ++it)
        if (it->second.members.size() > 0){
            S1 = it->second.members.size();
            H1 += (double)S1 * log((double)S1 / (double)S);
        }
            
    for (Groups::iterator it = (*gt).begin(); it != (*g).end(); ++it)
        if (it->second.members.size() > 0){
            S2 = it->second.members.size();
            H2 += (double)S2 * log((double)S2 / (double)S);
        }

    for (Groups::iterator it1 = (*g).begin(); it1 != (*g).end(); ++it1)
        for (Groups::iterator it2 = (*gt).begin(); it2 != (*gt).end(); ++it2)
            for(GroupNodes::iterator nit1 = it1->second.members.begin(); nit1 != it1->second.members.end(); ++nit1)
                for(GroupNodes::iterator nit2 = it2->second.members.begin(); nit2 != it2->second.members.end(); ++nit2)
                    if (nit1->second == nit2->second)
                        S12++;

    if (S12 > 0)
        H12 += (double)S12 * log((double)(S12 * S) / (double)(S1 * S2));

    return -2.0 * H12 / (H1 + H2);

}



/* ##################################### */
double group2GroupH(Group *g1, Group *g2, int K, GGLinks *gglinks){
	double H = 0.0;

	/* H += gsl_sf_lnfact(g1->g2glinks[g2->get_id()][0] + K -1); */
	H += gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][0] + K - 1);
	for(int i=1; i<K+1; ++i){
		/* H -= gsl_sf_lnfact(g1->g2glinks[g2->get_id()][i]); */
        H -= gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][i]);
	}
    return H;
}


/* ##################################### */
double hkState(int K, Groups *g1, Groups *g2, int d1_size, int d2_size, GGLinks *gglinks){
	double H = 0.0;

	for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
		for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
			H += group2GroupH(&(it1->second), &(it2->second), K, gglinks);
		
	H -= gsl_sf_lnfact(d1_size - g1->size()) + gsl_sf_lnfact(d2_size - g2->size());
	
	return H;

}


/* ##################################### */
void createRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, gsl_rng *rgen, int K, bool random, GGLinks *gglinks){
	int nweight[K+1];
	int ngrouplinks;
	int group;

	for (int i = 1; i<(*d1).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d1).size())) + 1) : group=i;
		(*d1)[i].setGroup(group);
		((*groups1)[group]) = Group(group, K);
		(*groups1)[group].members[i]=&((*d1)[i]);
	}
	
	for (int i = 1; i<(*d2).size()+1; ++i){
		(random) ? (group = (floor(gsl_rng_uniform(rgen) * (double)(*d2).size())) + 1) : group=i;
		(*d2)[i].setGroup(group);
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
		                if (nit->second.getId() == it2->second->getId()){
        		            ++ngrouplinks;
                            nweight[nit->second.getWeight()]++;
                        }
        		}

            (*gglinks)[g1->second.getId()][g2->second.getId()][0] = ngrouplinks;
		    for(int i=1; i<K+1; ++i)
                (*gglinks)[g1->second.getId()][g2->second.getId()][i] = nweight[i];
		}
	}
}

/* ##################################### */

double* genLogFactList(int size){
    double* logFactList = (double*) calloc(size, sizeof(double));

    for (int i = 0; i<size; i++)
        logFactList[i] = gsl_sf_lnfact(i);

    return logFactList;
}

double logFact(int key, int size, double* logFactList){
    if (size<key)
        return gsl_sf_lnfact(key);
    else
        return logFactList[key];
}

int mcStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinks){


    /* struct timeval stop, start; */
    /* gettimeofday(&start, NULL); */
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
    double set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);


    /* for (move=0; move<(nnod1+nnod2)*factor; move++) { */
        if (gsl_rng_uniform(stepgen) < set_ratio){
            g = g1;d_move = d1;d_nomove = d2;set_ind = true;
        }else{
            g = g2;d_move = d2;d_nomove = d1;set_ind = false;
        }
        set_size_move = d_move->size();


        dice = gsl_rng_uniform(stepgen) * set_size_move + 1;
        n = &(d_move->at(dice));
        oldgrp = n->getGroup();
        newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;
        while ( newgrp == oldgrp ) 
            newgrp = gsl_rng_uniform(groupgen) * set_size_move + 1;

        src_g = &(*g)[oldgrp];
        dest_g = &(*g)[newgrp];

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (!visitedgroup[id]){
                if (set_ind){
                    dH -= logFact((*gglinks)[src_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                    dH -= logFact((*gglinks)[dest_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                }else{
                    dH -= logFact((*gglinks)[id][src_g->getId()][0] + K - 1, logsize, lnfactlist);
                    dH -= logFact((*gglinks)[id][dest_g->getId()][0] + K - 1, logsize, lnfactlist);
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH -= -logFact((*gglinks)[src_g->getId()][id][i], logsize, lnfactlist);
                        dH -= -logFact((*gglinks)[dest_g->getId()][id][i], logsize, lnfactlist);
                    }else{
                        dH -= -logFact((*gglinks)[id][src_g->getId()][i], logsize, lnfactlist);
                        dH -= -logFact((*gglinks)[id][dest_g->getId()][i], logsize, lnfactlist);
                    }
                }
                visitedgroup[id] = true;
            }
        }

        if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
            dH += logFact(set_size_move - g1->size(), logsize, lnfactlist);

        if (set_ind){
            src_g->removeNodeS1(n, d_nomove, gglinks);
            dest_g->addNodeS1(n, d_nomove, gglinks);
        }else{
            src_g->removeNodeS2(n, d_nomove, gglinks);
            dest_g->addNodeS2(n, d_nomove, gglinks);
        }

        for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
            id = (*d_nomove)[it->second.getId()].getGroup();
            if (visitedgroup[id]){
                if (set_ind){
                    dH += logFact((*gglinks)[src_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                    dH += logFact((*gglinks)[dest_g->getId()][id][0] + K - 1, logsize, lnfactlist);
                }else{
                    dH += logFact((*gglinks)[id][src_g->getId()][0] + K - 1, logsize, lnfactlist);
                    dH += logFact((*gglinks)[id][dest_g->getId()][0] + K - 1, logsize, lnfactlist);
                }

                for (int i = 1; i<K+1; ++i){
                    if (set_ind){
                        dH += -logFact((*gglinks)[src_g->getId()][id][i], logsize, lnfactlist);
                        dH += -logFact((*gglinks)[dest_g->getId()][id][i], logsize, lnfactlist);
                    }else{
                        dH += -logFact((*gglinks)[id][src_g->getId()][i], logsize, lnfactlist);
                        dH += -logFact((*gglinks)[id][dest_g->getId()][i], logsize, lnfactlist);
                    }
                }
                visitedgroup[id] = false;
            }
        }	

        if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
            dH -= logFact(set_size_move - g1->size(), logsize, lnfactlist);

        if ( dH <= 0.0 || gsl_rng_uniform(stepgen) < exp(-dH)){
            *H += dH;
        }else{
            if(set_ind){
                dest_g->removeNodeS1(n, d_nomove, gglinks);
                src_g->addNodeS1(n, d_nomove, gglinks);
            }else{
                dest_g->removeNodeS2(n, d_nomove, gglinks);
                src_g->addNodeS2(n, d_nomove, gglinks);
            }
        }
        /* gettimeofday(&stop, NULL); */
        /* printf("Time %lu\n", stop.tv_usec - start.tv_usec); */
    /* } //End of MC step */


	return 0;
}

int getDecorrelationKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *stepgen, gsl_rng *groupgen, 
                double *H, int K, double *lnfactlist, int logsize, int nnod1, int nnod2, GGLinks *gglinksGroups){

    int nrep = 10, step, norm = 0;
    int x1, x2, i;
    double *decay, *decay1, *decay2;
    double meanDecay, meanDecay1, meanDecay2, sigmaDecay;
    double y11 = 0.0, y12 = 0.0, y21 = 0.0, y22 = 0.0, result;
	Groups g2t, g1t;

    decay1 = (double*) calloc(nrep, sizeof(double));
    decay2 = (double*) calloc(nrep, sizeof(double));

    x2 = (nnod1 + nnod2) / 5;
    if (x2 < 10)
        x2 = 10;
    x1 = x2 / 4;

    for (i=0; i<nrep; i++){

        printf("Estimating decorrelation time %d/%d\n", i, nrep);
        g1t = (*g1);
        g2t = (*g2);

        for (step=0; step<=x2; step++){
            mcStepKState(g1, g2, d1, d2, stepgen, groupgen, H, K, lnfactlist, logsize, nnod1, nnod2, gglinksGroups);

            if (step == x1){
                y11 = mutualInfo(g1, &g1t, nnod1, nnod2);
                y12 = mutualInfo(g2, &g2t, nnod1, nnod2);
            }
        }
        y21 = mutualInfo(g1, &g1t, nnod1, nnod2);
        y22 = mutualInfo(g2, &g2t, nnod1, nnod2);

        (nnod1>1) ? decay1[i] = 2 * getDecay(nnod1, x1, x2, y11, y21) : decay1[i] = 1.e-6;
        (nnod2>1) ? decay2[i] = 2 * getDecay(nnod2, x1, x2, y12, y22) : decay2[i] = 1.e-6;

        printf("# %d %d %d %d %f %f %f %f\n", nnod1, nnod2, x1, x2, y11, y21, y12, y22);
        printf("# Decorrelation times (estimate %d) = %g %g\n", i + 1, decay1[i], decay2[i]);

        if (decay1[i] < 0. || decay2[i] < 0.){
            i--;
            printf("\t# ignoring...\n");
        }
    }

    meanDecay1 = mean(decay1, nrep);
    meanDecay2 = mean(decay2, nrep);

    if (meanDecay1 > meanDecay2){
        meanDecay = meanDecay1;
        decay = decay1;
        sigmaDecay = stddev(decay1, nrep);
    }else{
        meanDecay = meanDecay2;
        decay = decay2;
        sigmaDecay = stddev(decay2, nrep);
    }

    result = meanDecay*nrep;

    for (i=0; i<nrep; i++)
       (fabs(decay[i] - meanDecay) / sigmaDecay > 2) ? result -= decay[i] : ++norm;

    printf("# Decorrelation step: %d\n", (int)(result / norm + 0.5));

    return (int)(result / norm + 0.5);
}


void printGroups(Groups g, int mark){
    for (Groups::iterator it = (g).begin(); it != (g).end(); ++it){
        std::cout << "[Group:" << it->second.getId();
        for(GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << ", Node:" << it1->second->getId() << " Links: ";
                for (Links::iterator nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); ++nit)
                    std::cout << nit->second.getId() << ", ";

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
    
    lnfactlist = genLogFactList(logsize);


	d1c = d1;
	d2c = d2;	
    nnod1 = d1.size();
    nnod2 = d2.size();

    GGLinks gglinks(boost::extents[nnod1+2][nnod2+2][mark+1]);



	createRandomGroups(&d1c, &d2c, &groups1, &groups2, groups_randomizer, mark, false, &gglinks);


    decorStep = getDecorrelationKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks);

    /* printf("DECORRELATION: %d\n", decorStep); */

	H = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks);
    double TH;
	std::cout << "Initial H: "<< H <<"\n";
	for(int i=0; i<STEPS; i++){
		mcStepKState(&groups1, &groups2, &d1c, &d2c, step_randomizer, groups_randomizer, &H, mark, lnfactlist, logsize, nnod1, nnod2, &gglinks);
        TH = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks);
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

