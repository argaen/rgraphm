#ifndef GROUP_H
#define GROUP_H

#include <vector>
#include "Node.h"
#include <boost/multi_array.hpp>

typedef boost::unordered_map<int, Node> Hash_Map;
typedef boost::unordered_map< int , Node* > GroupNodes;
typedef boost::multi_array<int, 3> GGLinks;

class Group {
	int id;
	int k;
	
	public:
		GroupNodes members;

		Group ();
		Group ( int id );
		Group ( int id, int K);

		int getId () const;
		void setId ( int id );

		int addNodeS1 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int addNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int removeNodeS1 ( Node *n, Hash_Map *d, GGLinks *gglinks );
		int removeNodeS2 ( Node *n, Hash_Map *d, GGLinks *gglinks );
};

#endif
