/* #include <iostream> */
#include <string.h>
#include <fstream>
/* #include <cstdlib> */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <boost/unordered_map.hpp>
/* #include <vector> */

#include "Node.h"
#include "Link.h"
#include "Group.h"

#define STEPS 100

typedef boost::unordered_map< int, Node > Nodeset;
typedef std::pair< int, Node > Nodeset_Pair;

typedef std::pair< int, Link > Link_Pair;

typedef std::pair< int, Node* > Member_Pair;

typedef boost::unordered_map< int, Group > Groups;


void parseArguments ( int argc, char **argv, char** inFile, char** qFile, int* stepseed, int* groupseed, int* mark ) {

	if ( argc != 11 ) {
		fprintf (stderr, "Usage: main_recommender -q queryFile -t trainFile -s stepseed -g groupsseed -m mark \n");
        exit(1);
	}

	int c;
	while ( (c = getopt (argc, argv, "t:q:s:g:m:")) != -1 )
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
				if ( optopt == 't' || optopt == 'q' || optopt == 'm' || optopt == 's' || optopt == 'g' )
					fprintf( stderr, "Option -%c requires an argument.\n", optopt );
				else if ( isprint ( optopt ) )
					fprintf( stderr, "Unknown option `-%c'.\n", optopt );
				else
					fprintf( stderr, "Unknown option character `\\x%x'.\n", optopt );
				exit(1);

			default:
				fprintf( stderr, "Usage: main_recommender -q queryFile -t trainFile -s seed -m mark \n" );
				exit(1);
		}
}


void showNodeSet ( Nodeset n ){
    std::cout << "--------------\n";
    for ( Nodeset::iterator it = n.begin(); it != n.end(); ++it ){
        std::cout << "Node " << it->second.get_id() << ": ";
        for ( Links::iterator nit = it->second.neighbours.begin(); nit != it->second.neighbours.end(); ++nit )
            std::cout << nit->second.get_id() << " " << nit->second.get_weight() << " | ";

        std::cout << "\n";
    }
    std::cout << "--------------\n";
}

void showGroupsInfo ( Groups group ){
    

	std::cout << " ################ G R O U P S ############### \n";
    for ( Groups::iterator it = ( group ).begin(); it != ( group ).end(); ++it){
		std::cout << "[Group:" << it->second.get_id();
        std::cout << " Members: "<< it->second.members.size() << "] =";
        for( GroupNodes::iterator it1 = it->second.members.begin(); it1 != it->second.members.end(); ++it1){
            std::cout << " Node:" << ( *( it1->second ) ).get_id() << " Links: ";
                for (Links::iterator nit = ( *( it1->second ) ).neighbours.begin(); nit != ( *( it1->second ) ).neighbours.end(); ++nit)
                    std::cout << nit->second.get_id() << ", ";
                std::cout <<" | ";

        }
		/* std::cout << "] \nTotalLinks: " <<it->second.get_nlinks() << "\n"; */
		/* std::cout << "Kratings: "; */
		/* /1* for ( int i=1; i<mark+1; ++i ) *1/ */
		/* /1* 	std::cout << it->second.kratings[i] << " "; *1/ */
		std::cout << "\n";
			
	}

}

int fillNodesets ( char* tFileName, Nodeset *nodeset1, Nodeset *nodeset2 ){

    int id1, id2, weight;

    std::ifstream inFile ( tFileName );
    if ( inFile.is_open() ){
        while (inFile >> id1 >> id2 >> weight){
            ( *nodeset1 ).insert( Nodeset_Pair( id1, Node( id1 ) ) );
            ( *nodeset1 )[id1].neighbours.insert( Link_Pair( id2, Link( id2, weight ) ) );

            ( *nodeset2 ).insert( Nodeset_Pair( id2, Node( id2 ) ) );
            ( *nodeset2 )[id2].neighbours.insert( Link_Pair( id1, Link( id1, weight ) ) );
        }
        inFile.close();
        return 0;
    }else{
        std::cout << "Couldn't open file " << tFileName << "\n";
        exit(1);
    }


}


void CreateRandomGroups ( Nodeset *nodeset1, Nodeset *nodeset2, Groups *groups1, Groups *groups2, gsl_rng *rgen, int K, bool random ){
    
    int nweight[K+1];
    int nlinks = 0, ngrouplinks = 0;
    int group;

    memset( nweight, 0, sizeof(int) * (K+1) );

    //Place nodes to groups randomly (or not)
    for ( int i = 1; i < ( *nodeset1 ).size() + 1; ++i ){
        ( random ) ? ( group = ( floor ( gsl_rng_uniform ( rgen ) * ( double )( *nodeset1 ).size() ) ) + 1 ) : group = i;
        ( *nodeset1 )[i].set_group( group );
        if ( ( ( *groups1 ).find(group) ) == ( *groups1 ).end() )
            ( *groups1 )[group] = Group( group, K );

        ( *groups1 )[group].members.insert( Member_Pair ( ( *nodeset1 )[i].get_id() , &( ( *nodeset1 )[i] )));

    }

    for ( int i = 1; i < ( *nodeset2 ).size() + 1; ++i ){
        ( random ) ? ( group = ( floor ( gsl_rng_uniform ( rgen ) * ( double )( *nodeset2 ).size() ) ) + 1 ) : group = i;
        ( *nodeset2 )[i].set_group( group );
        if ( ( ( *groups2 ).find(group) ) == ( *groups2 ).end() )
            ( *groups2 )[group] = Group( group, K );

        ( *groups2 )[group].members[( ( *nodeset2 )[i].get_id() )] = &( ( *nodeset2 )[i] );
    }

   
    //Fill Groups information (link counts, etc)

    for ( Groups::iterator it1 = groups1->begin(); it1 != groups1->end(); ++it1 ){
        for ( Groups::iterator it2 = groups2->begin(); it2 != groups2->end(); ++it2 ){
    
        }
    }



}


int main ( int argc, char **argv ){

    int mark;

    int groupseed, stepseed;
    gsl_rng *groups_randomizer, *step_randomizer;

    char *tFileName, *qFileName;

    Nodeset nodeset1, nodeset2;
    Groups groups1, groups2;

	parseArguments ( argc, argv, &tFileName, &qFileName, &stepseed, &groupseed, &mark );

    groups_randomizer = gsl_rng_alloc ( gsl_rng_mt19937 );
    step_randomizer = gsl_rng_alloc ( gsl_rng_mt19937 );
    gsl_rng_set ( step_randomizer, stepseed );
    gsl_rng_set ( groups_randomizer, groupseed );


    std::cout << stepseed << " " << tFileName << " " << qFileName << " " << groupseed << " " << mark << "\n";

    fillNodesets ( tFileName, &nodeset1, &nodeset2 );
    showNodeSet ( nodeset1 );
    showNodeSet ( nodeset2 );


    CreateRandomGroups ( &nodeset1, &nodeset2, &groups1, &groups2, groups_randomizer, mark, true );

    showGroupsInfo ( groups1 );
    showGroupsInfo ( groups2 );



}
