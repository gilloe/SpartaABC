#pragma once
#include<iostream>
//copied from my previous version OI 15.3

typedef unsigned short nucID;
class Node
{
private:
	nucID data;
	Node* next;
public:
	Node();
	Node(const nucID& data);
	Node(const nucID& data, Node& n);
	~Node();
	bool isEmpty();
	nucID getData();
	Node* getNext();
	void setNext(Node* n);
	void setData(const nucID& data);
};


