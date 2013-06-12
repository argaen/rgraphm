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


int Group::getId () const {

	return Group::id;
}


void Group::setId ( int id ) {

	Group::id = id;
}

int Group::addNodeS1 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	n->setGroup(Group::id);
	Group::members[n->getId()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        ++(*gglinks)[Group::id][d->at(it->second.getId()).getGroup()][0];
        ++(*gglinks)[Group::id][d->at(it->second.getId()).getGroup()][it->second.getWeight()];
	}
 	
	return 0;
}

int Group::addNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	n->setGroup(Group::id);
	Group::members[n->getId()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        ++(*gglinks)[d->at(it->second.getId()).getGroup()][Group::id][0];
        ++(*gglinks)[d->at(it->second.getId()).getGroup()][Group::id][it->second.getWeight()];
	}
	return 0;
}


int Group::removeNodeS1 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	Group::members.erase(n->getId());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        --(*gglinks)[Group::id][d->at(it->second.getId()).getGroup()][it->second.getWeight()];
        --(*gglinks)[Group::id][d->at(it->second.getId()).getGroup()][0];
	}
	return 0;
}

int Group::removeNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

	Group::members.erase(n->getId());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        --(*gglinks)[d->at(it->second.getId()).getGroup()][Group::id][it->second.getWeight()];
        --(*gglinks)[d->at(it->second.getId()).getGroup()][Group::id][0];
	}
	return 0;
}
