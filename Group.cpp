#include "Group.h"
#include <string.h>


Group::Group( int id ) {

	Group::id = id;
}

Group::Group ( int id, int K ) {
		
	Group::id = id;
	Group::k = K;

}

Group::Group() {}


int Group::get_id () const {

	return Group::id;
}


void Group::set_id ( int id ) {

	Group::id = id;
}

int Group::add_node_s1 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	n->set_group(Group::id);
	Group::members[n->get_id()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        ++(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][0];
        ++(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][it->second.get_weight()];
	}
 	
	return 0;
}

int Group::add_node_s2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	n->set_group(Group::id);
	Group::members[n->get_id()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        ++(*gglinks)[d->at(it->second.get_id()).get_group()][Group::id][0];
        ++(*gglinks)[d->at(it->second.get_id()).get_group()][Group::id][it->second.get_weight()];
	}
	return 0;
}


int Group::remove_node_s1 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	Group::members.erase(n->get_id());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        --(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][it->second.get_weight()];
        --(*gglinks)[Group::id][d->at(it->second.get_id()).get_group()][0];
	}
	return 0;
}

int Group::remove_node_s2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	Group::members.erase(n->get_id());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        --(*gglinks)[d->at(it->second.get_id()).get_group()][Group::id][it->second.get_weight()];
        --(*gglinks)[d->at(it->second.get_id()).get_group()][Group::id][0];
	}
	return 0;
}





