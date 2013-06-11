#ifndef GROUP_H
#define GROUP_H

#include <vector>
#include "Node.h"
#include <boost/multi_array.hpp>

typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map< int , Node* > GroupNodes;
typedef boost::multi_array<int, 3> GGLinks;
typedef GGLinks::index GG_index;

class Group {
	int id;
	int k;
	
	public:
		GroupNodes members;

		Group ();
		Group ( int id );
		Group ( int id, int K);

		int get_id () const;
		void set_id ( int id );

		int add_node_s1 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int add_node_s2 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int remove_node_s1 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int remove_node_s2 ( Node *n, Hash_Map *d, GGLinks *gglinks );
};

#endif
