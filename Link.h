#ifndef LINK_H
#define LINK_H

class Link {
	int id;
	int weight;

	public:
		Link();
		Link(int id, int weight);

		int get_id() const;
		int get_weight() const;

		void set_id(int id);
		void set_weight(int weigth);
};

#endif
