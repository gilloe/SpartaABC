#pragma once
#include<iostream>
#include "Node.h"
class LinkedList
{
private:
	Node* head;
	Node* tail;
	size_t length;
public:
	LinkedList();
	~LinkedList();

	void delete_LL();

	void print();
	size_t getLength();
	Node* getTail();
	Node* getHead();
	void add_list_to_list(LinkedList &ll, Node* insertafter);
	void addnode_to_list(Node* insertnode, Node* insertafter);
};

