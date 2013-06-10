#include "Group.h"
#include <string.h>


Group::Group( int id ) {

	Group::id = id;
}

Group::Group ( int id, int K ) {
		
	Group::id = id;
	Group::k = K;
    kratings = (int*) malloc((K+1)*sizeof(int));
    memset(kratings, 0, sizeof(int)*(K+1));

}

Group::Group() {}


int Group::get_id () const {

	return Group::id;
}


void Group::set_id ( int id ) {

	Group::id = id;
}


int Group::get_nlinks () const {
	
	return Group::nlinks;
}


void Group::set_nlinks ( int nlinks ) {

	Group::nlinks = nlinks;
}


int Group::add_node ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	n->set_group(Group::id);
	Group::members[n->get_id()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        ++(*gglinks)[Group::id][it->second.get_id()][it->second.get_weight()];
        ++(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][0];
	}
 	
	return 0;

}


int Group::remove_node ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	Group::members.erase(n->get_id());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        --(*gglinks)[Group::id][it->second.get_id()][it->second.get_weight()];
        --(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][0];
	}
	return 0;
}





