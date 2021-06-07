#include "Node.h"
//copied from my previous version OI 15.3
#pragma region ctors

Node::Node()
{
	data = NULL;
	next = NULL;
}

Node::Node(const nucID& data)
{
	this->data = data;
	next = NULL;
}

Node::Node(const nucID& data, Node& n)
{
	this->data = data;
	this->next = &n;
}
#pragma endregion
#pragma region dtors
Node::~Node()
{
}
#pragma endregion
#pragma region funcs
Node* Node::getNext()
{
	return next;
}

nucID Node::getData()
{
	return data;
}

void Node::setData(const nucID& data)
{
	this->data = data;
}

void Node::setNext(Node* n)
{
	this->next = n;
}
bool Node::isEmpty() {
	if (this->data== NULL) {
		return true;
	}
	return false;
}
#pragma endregion

