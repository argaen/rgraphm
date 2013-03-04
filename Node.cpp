#include "Node.h"

Node::Node (int id) { Node::id=id; }
Node::Node() {}


int Node::get_id() const {
	return Node::id;
}

int Node::get_group() const {
	return Node::group;
}


void Node::set_id(int id) {
	Node::id = id;
}

void Node::set_group(int group) {
	Node::group = group;
}
