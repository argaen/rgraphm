#ifndef GROUP_H
#define GROUP_H

#include <vector>
#include "Node.h"
#include <boost/multi_array.hpp>

typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map< int , Node* > GroupNodes;
typedef boost::unordered_map< int, boost::unordered_map < int, int > > G2GRelations;
typedef boost::multi_array<double, 3> GGLinks;
typedef GGLinks::index GG_index;

class Group {
	int id;
	int nlinks;
	int k;
	
	public:
		G2GRelations g2glinks;
		int *kratings;
		GroupNodes members;

		Group ();
		Group ( int id );
		Group ( int id, int K);

		int get_id () const;
		int get_nlinks () const;

		void set_id ( int id );
		void set_nlinks ( int nlinks );

		int add_node ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int remove_node ( Node *n, Hash_Map *d, GGLinks *gglinks );
};

#endif
