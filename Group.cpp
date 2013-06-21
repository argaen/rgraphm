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

    int gid;
	n->setGroup(Group::id);
	Group::members[n->getId()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        gid = d->at(it->second.getId()).getGroup();
        ++(*gglinks)[Group::id][gid][0];
        ++(*gglinks)[Group::id][gid][it->second.getWeight()+1];
	}
 	
	return 0;
}

int Group::addNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

    int gid;
	n->setGroup(Group::id);
	Group::members[n->getId()] = n;
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        gid = d->at(it->second.getId()).getGroup();
        ++(*gglinks)[gid][Group::id][0];
        ++(*gglinks)[gid][Group::id][it->second.getWeight()+1];
	}
	return 0;
}


int Group::removeNodeS1 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

    int gid;
	Group::members.erase(n->getId());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        gid = d->at(it->second.getId()).getGroup();
        --(*gglinks)[Group::id][gid][it->second.getWeight()+1];
        --(*gglinks)[Group::id][gid][0];
	}
	return 0;
}

int Group::removeNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks ) {

    int gid;
	Group::members.erase(n->getId());
	for(Links::iterator it = n->neighbours.begin(); it != n->neighbours.end(); ++it){
        gid = d->at(it->second.getId()).getGroup();
        --(*gglinks)[gid][Group::id][it->second.getWeight()+1];
        --(*gglinks)[gid][Group::id][0];
	}
	return 0;
}
