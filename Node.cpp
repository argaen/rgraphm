#include "Node.h"

Node::Node (int id) { Node::id=id; }
Node::Node() {}


int Node::getId() const {
	return Node::id;
}

int Node::getGroup() const {
	return Node::group;
}


void Node::setId(int id) {
	Node::id = id;
}

void Node::setGroup(int group) {
	Node::group = group;
}
