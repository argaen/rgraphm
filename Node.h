#ifndef NODE_H
#define NODE_H

#include "Link.h"
#include <boost/unordered_map.hpp>

typedef boost::unordered_map<int, Link> Links;

class Node {

	int id;
	int group;

	public:
		Links neighbours;
		Node();
		Node(int id);
		int get_id() const;
		int get_group() const;
		void set_id(int id);
		void set_group(int id);

};

#endif
