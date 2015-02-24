#include <iostream>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <vector>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <numeric>

#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>


#include "Node.h"
#include "Link.h"
#include "Group.h"
#include "utils.h"

#define LOGSIZE 5000

typedef boost::random::mt19937 RNGType;
RNGType rng;  //Constant output
/* RNGType rng(time(0)); */
typedef boost::unordered_map<double, double> HashExp;
typedef boost::unordered_map<int, double> LnFactList;
typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map<Node*, Node*> Queries;
typedef boost::unordered_map<std::string, int> IdMap;
typedef boost::multi_array<int, 3> GGLinks;
typedef boost::multi_array<double, 2> Scores;
typedef boost::numeric::ublas::vector<int> IVector;
typedef std::vector< boost::tuple<Node*, Node*> > tuple_list;



/*!
 * Calculate H between two Node groups. 
 * @param K Number of values a Link can be weighted.
 * @param g2glinks Pointer within the g2glinks matrix pointing to the specific values between this two Groups.
 * @return H between the two Groups.
 */
double group2GroupH(Group *g1, Group *g2, int K, GGLinks *gglinks) {
    double H = 0.0;

    /* H += gsl_sf_lnfact(g1->g2glinks[g2->get_id()][0] + K -1); */
    H += gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][0] + K - 1);
    for(int i=1; i<K+1; ++i) {
        /* H -= gsl_sf_lnfact(g1->g2glinks[g2->get_id()][i]); */
        H -= gsl_sf_lnfact((*gglinks)[g1->getId()][g2->getId()][i]);
    }
    return H;
}


/*
 * Calculate the H of all system at the current state. Consists on calculating the H between all the Group pairs (see group2groupH()). Used to get initial H of the system.
 * @param K Number of values a Link can be weighted.
 * @param g1 Set of Groups from the first partition.
 * @param g2 Set of Groups from the second partition.
 * @param d1_size Number of Nodes of the first partition.
 * @param d2_size Number of Nodes of the second partition.
 * @param g2glinks Matrix where all the Link count between all group pairs is stored.
 * @param ng1 Number of Groups of the first partition.
 * @param ng2 Number of Groups of the second partition.
 * @return Total H of the system.
 */
double hkState(int K, Groups *g1, Groups *g2, int d1_size, int d2_size, GGLinks *gglinks, int ng1, int ng2) {
    double H = 0.0;
    /* typedef GGLinks::index index; */

    for (Groups::iterator it1 = g1->begin(); it1 != g1->end(); ++it1)
        for (Groups::iterator it2 = g2->begin(); it2 != g2->end(); ++it2)
            H += group2GroupH(&(it1->second), &(it2->second), K, gglinks);

    H -= gsl_sf_lnfact(d1_size - ng1) + gsl_sf_lnfact(d2_size - ng2);

    /* for(index i = 0; i != 2; ++i){ */
    /*     for(index j = 0; j != 2; ++j){ */
    /*         std::cout << "|"; */
    /*         for(index k = 0; k != 3; ++k){ */
    /*             std::cout << (*gglinks)[i][j][k]; */
    /*         } */
    /*         std::cout << "   "; */
    /*     } */
    /*     std::cout << std::endl << std::endl; */
    /* } */

    return H;

}


/*
 * Create groups randomly from two different partitions.
 * @param d1 First partition of Nodes.
 * @param d2 Second partition of Nodes.
 * @param groups1 Structure where the first random set of Groups is stored.
 * @param groups2 Structure where the second random set of Groups is stored.
 * @K is the number of different weight values Links can have.
 * @gglinks 3D matrix where Link to Link weights are stored.
 * @nnod1 Total number of nodes in partition 1.
 * @nnod2 Total number of nodes in partition 2.
 */
void createRandomGroups(Hash_Map *d1, Hash_Map *d2, Groups *groups1, Groups *groups2, int K, GGLinks *gglinks, int nnod1, int nnod2) {
    int nweight[K];
    int ngrouplinks;
    int group, i = 0;
    Groups::iterator g1, g2;
    GroupNodes::iterator it1, it2;
    Links::iterator nit;
    Hash_Map::iterator itn;

    i = 0;
    for (itn = d1->begin(); itn != d1->end(); itn++) {
        group = i;
        itn->second.setGroup(group);
        ((*groups1)[group]) = Group(group, K);
        (*groups1)[group].members[itn->second.getId()]=&(itn->second);
        i++;
    }

    i = 0;
    for (itn = d2->begin(); itn != d2->end(); itn++) {
        group = i;
        itn->second.setGroup(group);
        ((*groups2)[group]) = Group(group, K);
        (*groups2)[group].members[itn->second.getId()]=&(itn->second);
        i++;
    }


    for (g1 = groups1->begin(); g1 != groups1->end(); g1++) {
        for (g2 = groups2->begin(); g2 != groups2->end(); g2++) {
            memset(nweight, 0, sizeof(int)*K);
            ngrouplinks = 0;

            //Groups 1 info filling

            for (it1 = g1->second.members.begin(); it1 != g1->second.members.end(); it1++)
                for (nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); nit++) {
                    for (it2 = g2->second.members.begin(); it2 != g2->second.members.end(); it2++)
                        if (nit->second.getId() == it2->second->getId()) {
                            ++ngrouplinks;
                            nweight[nit->second.getWeight()]++;
                        }
                }

            (*gglinks)[g1->second.getId()][g2->second.getId()][0] = ngrouplinks;
            for(i=1; i<K+1; ++i) {
                (*gglinks)[g1->second.getId()][g2->second.getId()][i] = nweight[i-1];
            }
        }
    }
}

double calculatedH(Node *n, Hash_Map *d_move, Hash_Map *d_nomove, GGLinks *gglinks, Group *src_g, Group *dest_g, int K, double *lnfactlist, bool set_ind, int *ng) {
    double dH = 0.0;
    int id;
    bool visitedgroup[d_nomove->size()];
    int set_size_move = d_move->size();

    memset(visitedgroup, false, d_nomove->size());

    for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it) {
        id = (*d_nomove)[it->second.getId()].getGroup();
        if (!visitedgroup[id]) {
            if (set_ind) {
                dH -= lnfactlist[(*gglinks)[src_g->getId()][id][0] + K - 1];
                dH -= lnfactlist[(*gglinks)[dest_g->getId()][id][0] + K - 1];
            } else {
                dH -= lnfactlist[(*gglinks)[id][src_g->getId()][0] + K - 1];
                dH -= lnfactlist[(*gglinks)[id][dest_g->getId()][0] + K - 1];
            }

            for (int i = 1; i<K+1; ++i) {
                if (set_ind) {
                    dH -= -lnfactlist[(*gglinks)[src_g->getId()][id][i]];
                    dH -= -lnfactlist[(*gglinks)[dest_g->getId()][id][i]];
                } else {
                    dH -= -lnfactlist[(*gglinks)[id][src_g->getId()][i]];
                    dH -= -lnfactlist[(*gglinks)[id][dest_g->getId()][i]];
                }
            }
            visitedgroup[id] = true;
        }
    }

    if ( src_g->members.size() == 1 || dest_g->members.size() == 0)
        dH += lnfactlist[set_size_move - *ng];

    if (set_ind) {
        src_g->removeNodeS1(n, d_nomove, gglinks);
        dest_g->addNodeS1(n, d_nomove, gglinks);
    } else {
        src_g->removeNodeS2(n, d_nomove, gglinks);
        dest_g->addNodeS2(n, d_nomove, gglinks);
    }

    if (src_g->members.size() == 0) *ng -=1;
    if (dest_g->members.size() == 1) *ng +=1;

    // Calculations for checking the dH when adding the node
    for (Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it) {
        id = (*d_nomove)[it->second.getId()].getGroup();
        if (visitedgroup[id]) {
            if (set_ind) {
                dH += lnfactlist[(*gglinks)[src_g->getId()][id][0] + K - 1];
                dH += lnfactlist[(*gglinks)[dest_g->getId()][id][0] + K - 1];
            } else {
                dH += lnfactlist[(*gglinks)[id][src_g->getId()][0] + K - 1];
                dH += lnfactlist[(*gglinks)[id][dest_g->getId()][0] + K - 1];
            }

            for (int i = 1; i<K+1; ++i) {
                if (set_ind) {
                    dH += -lnfactlist[(*gglinks)[src_g->getId()][id][i]];
                    dH += -lnfactlist[(*gglinks)[dest_g->getId()][id][i]];
                } else {
                    dH += -lnfactlist[(*gglinks)[id][src_g->getId()][i]];
                    dH += -lnfactlist[(*gglinks)[id][dest_g->getId()][i]];
                }
            }
            visitedgroup[id] = false;
        }
    }

    if ( src_g->members.size() == 0 || dest_g->members.size() == 1)
        dH -= lnfactlist[set_size_move - *ng];

    return dH;
}

double fastexp(double a, HashExp explist) {
    if (explist.find(a) == explist.end())
        explist[a] = exp(a);

    return explist[a];
}


/*!
 * Perform a Gibbs step in our system. It consists on looping over the nodes of both partitions calculating the dH for each possible movement. Once each dH is calculated for
 * a node, roll a dice and pick a movement and apply it.
 * @param g1 Set of Groups from the first partition.
 * @param g2 Set of Groups from the second partition.
 * @param d1 Set of Nodes from the first partition.
 * @param d2 Set of Nodes from the second partition.
 * @param rgen Random generator.
 * @param H Energy of the system.
 * @param lnfactlist Static list of factorial logarithm.
 * @param gglinks MAtrix where all the link count between all Group pairs is stored.
 * @param keys1 List of Node ids from the first partition.
 * @param keys2 List of Node ids from the second partition.
 * @param ng1 Number of (non empty) Groups from the first partition (g1).
 * @param ng2 Number of (non empty) Groups from the second partition (g2).
 */
int gibbsStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int nnod1, int nnod2, GGLinks *gglinks, int decorStep, IVector *keys1, IVector *keys2, int *ng1, int *ng2) {

    std::vector<std::tuple<int,double>> dh_values;
    std::vector<double> weights;
    Group *src_g, *dest_g;
    int newgrp, oldgrp, id_dest;
    Node *n;
    double dH = 0.0, empty_dH = 0.0;
    Groups tmp_groups;
    int empty_groups;
    HashExp explist;

    for (Hash_Map::iterator nodes1 = d1->begin(); nodes1 != d1->end(); nodes1++) {

        weights = {};
        dh_values = {};
        n = &(d1->at(nodes1->second.getId()));
        oldgrp = n->getGroup();
        src_g = &(*g1)[oldgrp];
        empty_groups = empty_dH = 0;

        for (Groups::iterator groups1 = g1->begin(); groups1 != g1->end(); groups1++) {
            newgrp = groups1->second.getId();
            if (newgrp == oldgrp) {
                dH = 0.0;
                dh_values.push_back(std::make_tuple(newgrp, 0.0));
                weights.push_back(1.0);

            } else {
                dest_g = &(*g1)[newgrp];
                if (dest_g->members.size() != 0 || empty_groups == 0) {
                    dH = calculatedH(n, d1, d2, gglinks, src_g, dest_g, K, lnfactlist, true, ng1);
                    dest_g->removeNodeS1(n, d2, gglinks);
                    src_g->addNodeS1(n, d2, gglinks);
                    if (src_g->members.size() == 1) *ng1 +=1;
                    if (dest_g->members.size() == 0) *ng1 -=1;
                    dh_values.push_back(std::make_tuple(newgrp, dH));
                    weights.push_back(fastexp(-dH, explist));

                    if (empty_groups == 0 && dest_g->members.size() == 0) {
                        empty_dH = dH;
                        empty_groups++;
                    }
                } else {
                    dh_values.push_back(std::make_tuple(newgrp, empty_dH));
                    weights.push_back(fastexp(-empty_dH, explist));
                }
            }
        }
        boost::random::discrete_distribution<int,double> weightdist(weights);
        boost::variate_generator<RNGType&, boost::random::discrete_distribution<int,double> > weightsampler(rng, weightdist);  //Constant output
        id_dest = weightsampler();
        dest_g = &((*g1)[std::get<0>(dh_values[id_dest])]);
        if (dest_g->getId() != src_g->getId()) {
            src_g->removeNodeS1(n, d2, gglinks);
            dest_g->addNodeS1(n, d2, gglinks);
            if (src_g->members.size() == 0) *ng1 -=1;
            if (dest_g->members.size() == 1) *ng1 +=1;
        }
        *H += std::get<1>(dh_values[id_dest]);
    }

    for (Hash_Map::iterator nodes2 = d2->begin(); nodes2 != d2->end(); nodes2++) {

        weights = {};
        dh_values = {};
        n = &(d2->at(nodes2->second.getId()));
        oldgrp = n->getGroup();
        src_g = &(*g2)[oldgrp];
        empty_groups = empty_dH = 0;

        for (Groups::iterator groups2 = g2->begin(); groups2 != g2->end(); groups2++) {
            newgrp = groups2->second.getId();
            if (newgrp == oldgrp) {
                dH = 0.0;
                dh_values.push_back(std::make_tuple(newgrp, 0.0));
                weights.push_back(1.0);

            } else {
                dest_g = &(*g2)[newgrp];
                if (dest_g->members.size() != 0 || empty_groups == 0) {
                    dH = calculatedH(n, d2, d1, gglinks, src_g, dest_g, K, lnfactlist, false, ng2);

                    dest_g->removeNodeS2(n, d1, gglinks);
                    src_g->addNodeS2(n, d1, gglinks);
                    if (src_g->members.size() == 1) *ng2 +=1;
                    if (dest_g->members.size() == 0) *ng2 -=1;
                    dh_values.push_back(std::make_tuple(groups2->second.getId(), dH));
                    weights.push_back(fastexp(-dH, explist));

                    if (empty_groups == 0 && dest_g->members.size() == 0) {
                        empty_dH = dH;
                        empty_groups++;
                    }
                } else {
                    dh_values.push_back(std::make_tuple(newgrp, empty_dH));
                    weights.push_back(fastexp(-empty_dH, explist));
                }

            }
        }
        boost::random::discrete_distribution<int,double> weightdist(weights);
        boost::variate_generator<RNGType&, boost::random::discrete_distribution<int,double> > weightsampler(rng, weightdist);  //Constant output
        id_dest = weightsampler();
        dest_g = &((*g2)[std::get<0>(dh_values[id_dest])]);
        if (dest_g->getId() != src_g->getId()) {
            src_g->removeNodeS2(n, d1, gglinks);
            dest_g->addNodeS2(n, d1, gglinks);
            if (src_g->members.size() == 0) *ng2 -=1;
            if (dest_g->members.size() == 1) *ng2 +=1;
        }
        *H += std::get<1>(dh_values[id_dest]);
    }

    return 0;
}


/*!
 * Perform a Monte Carlo step in our system. It consists on moving all the nodes from one or another partition (randomly chosen) decorrelationStep times between groups of the
 * same partition. Each time we calculate dH and if its valid (we want to minimize H), we update H value and continue with next step. If not, we undo the movement and continue with
 * the next step.
 * @param g1 Set of Groups from the first partition.
 * @param g2 Set of Groups from the second partition.
 * @param d1 Set of Nodes from the first partition.
 * @param d2 Set of Nodes from the second partition.
 * @param rgen Random generator.
 * @param H Energy of the system.
 * @param lnfactlist Static list of factorial logarithm.
 * @param nnod1 Number of Nodes in first partition (d1).
 * @param nnod2 Number of Nodes in second partition (d2).
 * @param gglinks Matrix where all the link count between all Group pairs is stored.
 * @param decorStep Coeficient to make this Step significative for the system.
 * @param keys1 List of Node ids from the first partition.
 * @param keys2 List of Node ids from the second partition.
 * @param ng1 Number of (non empty) Groups from the first partition (g1).
 * @param ng2 Number of (non empty) Groups from the second partition (g2).
 */
int mcStepKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int nnod1, int nnod2, GGLinks *gglinks, int decorStep, IVector *keys1, IVector *keys2, int *ng1, int *ng2) {


    Groups *g;
    Hash_Map *d_move;    // Partition where the movement occurs
    Hash_Map *d_nomove;  // Partition where the movement does not occur

    int factor;
    Group *src_g, *dest_g;
    int newgrp, oldgrp, dice, set_size_move;
    bool set_ind;
    Node *n;
    double dH;
    int *ng;
    double set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);
    IVector *keys;


    factor = (nnod1+nnod2)*decorStep;
    for (int move=0; move<factor; move++) {
        dH = 0.0;
        if (gsl_rng_uniform(rgen) < set_ratio) {
            g = g1;d_move = d1;d_nomove = d2;set_ind = true;keys = keys1; ng = ng1;
        } else {
            g = g2;d_move = d2;d_nomove = d1;set_ind = false;keys = keys2; ng = ng2;
        }
        set_size_move = d_move->size();

        dice = gsl_rng_uniform(rgen) * set_size_move;
        n = &(d_move->at((*keys)[dice]));
        oldgrp = n->getGroup();
        newgrp = gsl_rng_uniform(rgen) * set_size_move;
        while ( newgrp == oldgrp )
            newgrp = gsl_rng_uniform(rgen) * set_size_move;

        src_g = &(*g)[oldgrp];
        dest_g = &(*g)[newgrp];

        // Calculations for checking the dH when removing the node
        dH = calculatedH(n, d_move, d_nomove, gglinks, src_g, dest_g, K, lnfactlist, set_ind, ng);

        if ( dH <= 0.0 || gsl_rng_uniform(rgen) < exp(-dH)) {
            *H += dH;
        } else {
            if(set_ind) {
                dest_g->removeNodeS1(n, d_nomove, gglinks);
                src_g->addNodeS1(n, d_nomove, gglinks);
            } else {
                dest_g->removeNodeS2(n, d_nomove, gglinks);
                src_g->addNodeS2(n, d_nomove, gglinks);
            }
            if (src_g->members.size() == 1) *ng +=1;
            if (dest_g->members.size() == 0) *ng -=1;
        }
    }

    return 0;
}


/*!
 * Equilibrate the system before we start the sampling.
 * @param g1 Set of Groups from the first partition.
 * @param g2 Set of Groups from the second partition.
 * @param d1 Set of Nodes from the first partition.
 * @param d2 Set of Nodes from the second partition.
 * @param rgen Random generator.
 * @param H Energy of the system.
 * @param lnfactlist Static list of factorial logarithm.
 * @param nnod1 Number of Nodes in first partition (d1).
 * @param nnod2 Number of Nodes in second partition (d2).
 * @param gglinks Matrix where all the Link count between all Group pairs is stored.
 * @param decorStep Coeficient to make this Step significative for the system.
 * @param keys1 List of Node ids from the first partition.
 * @param keys2 List of Node ids from the second partition.
 * @param ng1 Number of (non empty) Groups from the first partition (g1).
 * @param ng2 Number of (non empty) Groups from the second partition (g2).
 */
int thermalizeMCKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int nnod1, int nnod2, GGLinks *gglinks, int decorStep, IVector *keys1, IVector *keys2, int *ng1, int *ng2, int (*algorithm)(Groups*, Groups*, Hash_Map*, Hash_Map*, gsl_rng*, double*, int, double*, int, int, GGLinks*, int, IVector*, IVector*, int*, int*)) {

    double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
    int nrep=20;
    int equilibrated=0;

    Hvalues = (double*)calloc(nrep, sizeof(double));

    do {
        for (int i=0; i<nrep; i++) {
            algorithm(g1, g2, d1, d2, rgen, H, K, lnfactlist, nnod1, nnod2, gglinks, decorStep, keys1, keys2, ng1, ng2);
            Hvalues[i] = *H;
        }

        HMean1 = mean(Hvalues, nrep);
        HStd1 = stddev(Hvalues, nrep);

        if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) {
            equilibrated++;
            fprintf(stderr, "#\tequilibrated (%d/5) H=%lf\n", equilibrated, HMean1);

        } else {
            fprintf(stderr, "#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n", HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep));
            HMean0 = HMean1;
            HStd0 = HStd1;
            equilibrated = 0;
        }
    } while (equilibrated < 5);

    return 0;
}


/*!
 * Look for how many movements we have to do before dH is significative.
 * @param g1 Set of Groups from the first partition.
 * @param g2 Set of Groups from the second partition.
 * @param d1 Set of Nodes from the first partition.
 * @param d2 Set of Nodes from the second partition.
 * @param rgen Random generator.
 * @param H Energy of the system.
 * @param lnfactlist Static list of factorial logarithm.
 * @param nnod1 Number of Nodes in first partition (d1).
 * @param nnod2 Number of Nodes in second partition (d2).
 * @param gglinks Matrix where all the Link count between all Group pairs is stored.
 * @param keys1 List of Node ids from the first partition.
 * @param keys2 List of Node ids from the second partition.
 * @param ng1 Number of (non empty) Groups from the first partition (g1).
 * @param ng2 Number of (non empty) Groups from the second partition (g2).
 * @return decorrelation step.
 */
int getDecorrelationKState(Groups *g1, Groups *g2, Hash_Map *d1, Hash_Map *d2, gsl_rng *rgen, double *H, int K, double *lnfactlist, int nnod1, int nnod2, GGLinks *gglinks, IVector *keys1, IVector *keys2, int *ng1, int *ng2) {

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

    for (i=0; i<nrep; i++) {

        fprintf(stderr, "Estimating decorrelation time %d/%d\n", i+1, nrep);
        g1t = (*g1);
        g2t = (*g2);

        for (step=0; step<=x2; step++) {
            mcStepKState(g1, g2, d1, d2, rgen, H, K, lnfactlist, nnod1, nnod2, gglinks, 1, keys1, keys2, ng1, ng2);

            if (step == x1) {
                y11 = mutualInfo(g1, &g1t, nnod1, nnod2);
                y12 = mutualInfo(g2, &g2t, nnod1, nnod2);
            }
        }
        y21 = mutualInfo(g1, &g1t, nnod1, nnod2);
        y22 = mutualInfo(g2, &g2t, nnod1, nnod2);

        (nnod1>1) ? decay1[i] = 2 * getDecay(nnod1, x1, x2, y11, y21) : decay1[i] = 1.e-6;
        (nnod2>1) ? decay2[i] = 2 * getDecay(nnod2, x1, x2, y12, y22) : decay2[i] = 1.e-6;

        fprintf(stderr, "# Decorrelation times (estimate %d) = %g %g\n\n", i + 1, decay1[i], decay2[i]);

        if (decay1[i] < 0. || decay2[i] < 0.) {
            i--;
            fprintf(stderr, "\t# ignoring...\n");
        }
    }

    meanDecay1 = mean(decay1, nrep);
    meanDecay2 = mean(decay2, nrep);

    if (meanDecay1 > meanDecay2) {
        meanDecay = meanDecay1;
        decay = decay1;
        sigmaDecay = stddev(decay1, nrep);
    } else {
        meanDecay = meanDecay2;
        decay = decay2;
        sigmaDecay = stddev(decay2, nrep);
    }

    result = meanDecay*nrep;

    for (i=0; i<nrep; i++)
        (fabs(decay[i] - meanDecay) / sigmaDecay > 2) ? result -= decay[i] : ++norm;

    return (int)(result / norm + 0.5);
}




int main(int argc, char **argv) {

    int decorStep = 0;
    int weight, logsize=LOGSIZE;
    std::string id1, id2;
    double H;
    double *lnfactlist;
    tuple_list queries;
    Hash_Map d1, d2;
    IdMap rd1, rd2;
    Groups groups2, groups1;
    gsl_rng *randomizer;
    char* tFileName;
    char* qFileName;
    char* str_algorithm = (char*) malloc(25*sizeof(char));
    int stepseed;
    int iterations = 10000;
    int mark, k, i, j, n, nk, tmpid1 = 0, tmpid2 = 0, tid1, tid2;
    int nnod1, nnod2, nqueries = 0;
    int ng1 = 0, ng2 = 0;
    time_t start;
    time_t end;
    double TH;
    std::ofstream outfile;
    int (*algorithm)(Groups*, Groups*, Hash_Map*, Hash_Map*, gsl_rng*, double*, int, double*, int, int, GGLinks*, int, IVector*, IVector*, int*, int*) = mcStepKState;

    memcpy(str_algorithm, "metropolis", sizeof("metropolis"));
    parseArguments(argc, argv, &tFileName, &qFileName, &stepseed, &mark, &iterations, &str_algorithm);
    if (strcmp(str_algorithm, "gibbs") == 0)
        algorithm = gibbsStepKState;
    else if (strcmp(str_algorithm, "metropolis") == 0)
        algorithm = mcStepKState;



    randomizer = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(randomizer, stepseed);

    fprintf(stderr, "Reading input files and initializing data structures...\n");
    std::ifstream tFile(tFileName);
    if (tFile.is_open()) {
        while (tFile >> id1 >> id2 >> weight) {
            if (rd1.find(id1) == rd1.end()) {
                rd1.insert(IdMap::value_type(id1,tmpid1));
                tmpid1++;
            }
            if (rd2.find(id2) == rd2.end()) {
                rd2.insert(IdMap::value_type(id2,tmpid2));
                tmpid2++;
            }
            tid1 = rd1[id1]; 
            tid2 = rd2[id2];
            d1.insert(Hash_Map::value_type(tid1,Node(id1, tid1)));
            d1[tid1].neighbours.insert(Links::value_type(tid2,Link(tid2, weight)));
            d2.insert(Hash_Map::value_type(tid2,Node(id2, tid2)));
            d2[tid2].neighbours.insert(Links::value_type(tid1,Link(tid1, weight)));
        }
        tFile.close();
    } else
        std::cout << "Couldn't open file " << tFileName << "\n";

    std::ifstream qFile(qFileName);
    if (qFile.is_open()) {
        while (qFile >> id1 >> id2) {
            if (rd1.find(id1) == rd1.end()) {
                rd1.insert(IdMap::value_type(id1,tmpid1));
                tid1 = rd1[id1]; 
                d1.insert(Hash_Map::value_type(tid1,Node(id1, tid1)));
                tmpid1++;
            }

            if (rd2.find(id2) == rd2.end()) {
                rd2.insert(IdMap::value_type(id2,tmpid2));
                tid2 = rd2[id2];
                d2.insert(Hash_Map::value_type(tid2,Node(id2, tid2)));
                tmpid2++;
            }

            tid1 = rd1[id1]; 
            tid2 = rd2[id2];
            /* if (d1[tid1].neighbours.find(tid2) != d1[tid1].neighbours.end()) { */
            d1[tid1].neighbours.erase(tid2);
            d2[tid2].neighbours.erase(tid1);


            queries.push_back(boost::tuple<Node*, Node*>(&(d1[rd1[id1]]), &(d2[rd2[id2]])));
        }
        qFile.close();
    } else
        std::cout << "Couldn't open file " << tFileName << "\n";

    lnfactlist = genLogFactList(logsize);

    nnod1 = d1.size();
    nnod2 = d2.size();
    ng1 = d1.size();
    ng2 = d2.size();
    nqueries = queries.size();
    fprintf(stdout, "Network size:\n \tnodes1: %d\n\tnodes2: %d\n\tnqueries: %d\n", nnod1, nnod2, nqueries);

    IVector keys1(nnod1), keys2(nnod2);
    i = 0;
    for(Hash_Map::iterator it = d1.begin(); it != d1.end(); ++it) {
        keys1[i] = it->first;
        ++i;
    }
    i = 0;
    for(Hash_Map::iterator it = d2.begin(); it != d2.end(); ++it) {
        keys2[i] = it->first;
        ++i;
    }

    GGLinks gglinks(boost::extents[nnod1][nnod2][mark+1]);
    Scores scores(boost::extents[nqueries][mark]);

    fprintf(stderr, "Creating groups...\n\n");
    createRandomGroups(&d1, &d2, &groups1, &groups2, mark, &gglinks, nnod1, nnod2);
    /* printf("############GROUPS1################\n"); */
    /* printGroups(groups1, mark); */
    /* printf("############GROUPS2################\n"); */
    /* printGroups(groups2, mark); */
    /* printf("###################################\n"); */


    fprintf(stderr, "Calculating initial H...\n");
    start = time(NULL);
    H = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks, ng1, ng2);
    end = time(NULL);
    fprintf(stderr, "Initial H: %f | Time Spent: %lu secs\n\n", H, end-start);

    if (strcmp(str_algorithm, "metropolis") == 0) {
        fprintf(stderr, "Calculating decorrelation step...\n");
        start = time(NULL);
        decorStep = getDecorrelationKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, nnod1, nnod2, &gglinks, &keys1, &keys2, &ng1, &ng2);
        end = time(NULL);
        fprintf(stderr, "Decorrelation step: %d | Time Spent: %lu secs\n\n", decorStep, end-start);
    }

    fprintf(stderr, "Thermalizing...\n");
    start = time(NULL);
    thermalizeMCKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, nnod1, nnod2, &gglinks, decorStep, &keys1, &keys2, &ng1, &ng2, algorithm);
    end = time(NULL);
    fprintf(stderr, "Time Spent: %lu secs\n\n", end-start);



    fprintf(stderr, "Starting %s Steps...\n", str_algorithm);
    start = time(NULL);
    int indquery;
    for(i=0; i<iterations; i++) {
        algorithm(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, nnod1, nnod2, &gglinks, decorStep, &keys1, &keys2, &ng1, &ng2);
        /* gibbsStepKState(&groups1, &groups2, &d1, &d2, randomizer, &H, mark, lnfactlist, nnod1, nnod2, &gglinks, decorStep, &keys1, &keys2, &ng1, &ng2); */
        /* printf("############GROUPS1################\n"); */
        /* printGroups(groups1, mark); */
        /* printf("############GROUPS2################\n"); */
        /* printGroups(groups2, mark); */
        /* printf("###################################\n"); */

        std::cout << i << " " << H << "\n";
        /* TH = hkState(mark, &groups1, &groups2, nnod1, nnod2, &gglinks, ng1, ng2); */
        /* std::cout << std::setprecision(20) << H << "    " << TH << "\n\n"; */
        indquery = 0;
        for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
            for (k=1; k<mark+1; k++) {
                nk = gglinks[(*it->get<0>()).getGroup()][(*it->get<1>()).getGroup()][k];
                n = gglinks[(*it->get<0>()).getGroup()][(*it->get<1>()).getGroup()][0];
                scores[indquery][k-1] += (float)(nk + 1) / (float)(n + mark);
            }
            indquery++;
        }


        if (i % 100 == 0) { 
            outfile.open("scores.tmp");
            indquery = 0;
            for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
                outfile << "| " << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << " | ";
                for (k=0; k<mark; k++)
                    outfile << scores[indquery][k]/(double)iterations << " ";
                indquery++;
                outfile << "\n";
            }
            outfile.close();
        }
    }
    end = time(NULL);
    fprintf(stderr, "Done. Time Spent: %lu secs\n\n", end-start);

    for (j=0; j<nqueries; j++)
        for (k=0; k<mark; k++)
            scores[j][k] /= (double)iterations;

    outfile.open("scores.tmp");
    fprintf(stdout, "RESULTS:\n");
    indquery = 0;
    for (tuple_list::iterator it = queries.begin(); it != queries.end(); ++it) {
        outfile << "# " << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << " # ";
        std::cout << (*it->get<0>()).getRealId() <<' '<< (*it->get<1>()).getRealId() << '\n';
        for (k=0; k<mark; k++) {
            outfile << scores[indquery][k] << " ";
            std::cout << scores[indquery][k] << " ";
        }
        indquery++;
        outfile << "\n\n";
        std::cout << "\n\n";
    }

    free(lnfactlist);
    gsl_rng_free(randomizer);
    return 0;

}
