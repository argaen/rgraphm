#include "Node.h"

Node::Node (std::string id, int tmpid) { Node::realid=id; Node::id=tmpid; }
Node::Node() {}


int Node::getId() const {
	return Node::id;
}

std::string Node::getRealId() const {
	return Node::realid;
}

int Node::getGroup() const {
	return Node::group;
}


void Node::setId(int id) {
	Node::id = id;
}

void Node::setRealId(std::string id) {
	Node::realid = id;
}

void Node::setGroup(int group) {
	Node::group = group;
}
