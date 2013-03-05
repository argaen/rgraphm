#include "Group.h"
#include <string.h>


Group::Group( int id ) {

	Group::id = id;
}

Group::Group ( int id, int K ) {
		
	Group::id = id;
	Group::k = K;
	kratings = new int[K+1];
	memset( kratings, 0, (K+1)*sizeof(int) );

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


int Group::add_node ( Node n ) {

	/* n.set_group(Group::id); */
	
	/* Group::members[n.get_id()] = n; */
	/* Group::set_nlinks(Group::nlinks + n.neighbours.size()); */
	/* for(Links::iterator it = n.neighbours.begin(); it != n.neighbours.end(); ++it){ */
	/* 	++Group::g2glinks[it->second.get_id()][it->second.get_weight()]; */
	/* 	++Group::g2glinks[it->second.get_id()][0]; */
	/* //	++Group::kratings[it->second.get_weight()]; */
	/* } */
 	
	return 0;

}


int Group::remove_node ( Node n ) {

	/* Group::members.erase(n.get_id()); */
	/* Group::set_nlinks(Group::nlinks - n.neighbours.size()); */
	/* for(Links::iterator it = n.neighbours.begin(); it != n.neighbours.end(); ++it){ */
	/* 	--Group::g2glinks[it->second.get_id()][it->second.get_weight()]; */
	/* 	--Group::g2glinks[it->second.get_id()][0]; */
	/* 	//--Group::kratings[it->second.get_weight()]; */
	/* } */


	return 0;

}





