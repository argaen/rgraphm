#ifndef NODE_H
#define NODE_H

#include "Link.h"
#include "string.h"
#include <boost/unordered_map.hpp>

typedef boost::unordered_map<int, Link> Links;

class Node {
    int id;
    std::string realid;
	int group;

	public:
		Links neighbours;

		Node();
		Node(std::string id, int tmpid);

		int getId() const;
        std::string getRealId() const;
		int getGroup() const;

		void setId(int id);
		void setRealId(std::string id);
		void setGroup(int id);
};

#endif
