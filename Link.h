#ifndef LINK_H
#define LINK_H

class Link {
	int id;
	int weight;

	public:
		Link();
		Link(int id, int weight);

		int getId() const;
		int getWeight() const;

		void setId(int id);
		void setWeight(int weigth);
};

#endif
