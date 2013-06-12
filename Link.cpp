#include "Link.h"

Link::Link (int id, int weight) { Link::id=id;Link::weight=weight; }
Link::Link() {}

int Link::getId() const {
	return Link::id;
}

int Link::getWeight() const {
	return Link::weight;
}

void Link::setId(int id) {
	Link::id = id;
}

void Link::setWeight(int weigth) {
	Link::weight = weigth;
}
