#include "utils.h"
#include <iostream>

void parseArguments(int argc, char **argv, char** inFile, char** qFile, int* stepseed, int* mark) {
	if (argc != 9) {

		fprintf (stderr, "Usage: main_recommender -q queryFile -t trainFile -s randomseed -m mark \n");
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

				fprintf (stderr, "Usage: main_recommender -q queryFile -t trainFile -s randomseed -m mark \n");
				exit(1);
		}
}

int ExponentialRootF(const gsl_vector *params, void *points, gsl_vector *f) {
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

void printGroups(Groups g, int mark){
    for (Groups::iterator it = (g).begin(); it != (g).end(); ++it){
        /* std::cout << "[Group:" << it->second.getId(); */
        for(GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << it1->second->getId() << " " ;
                /* for (Links::iterator nit = it1->second->neighbours.begin(); nit != it1->second->neighbours.end(); ++nit) */
                /*     std::cout << nit->second.getId() << ", "; */

        }
        std::cout << "\n";
    }
}
