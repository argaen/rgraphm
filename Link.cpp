#include "Link.h"

Link::Link (int id, int weight) { Link::id=id;Link::weight=weight; }
Link::Link() {}

int Link::get_id() const {
	return Link::id;
}

int Link::get_weight() const {
	return Link::weight;
}

void Link::set_id(int id) {
	Link::id = id;
}

void Link::set_weight(int weigth) {
	Link::weight = weigth;
}
