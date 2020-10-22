#include "Block.h"
#pragma region Ctrs
Block::Block(int ind, int len) {
	this->index = ind;
	this->length = len;
}
Block::Block()
{
}
#pragma endregion
#pragma region Dtrs

Block::~Block()
{
}
#pragma endregion
#pragma region funcs
bool Block::operator==(Block block)
{
	if (&block == this)
	{
		return true;
	}
	return false;
}

int Block::getIndex() {
	return this->index;
}
int Block::getLength() {
	return this->length;
}
void Block::setIndex(int ind) {
	this->index = ind;
}
void Block::setLength(int len) {
	this->length = len;
}


#pragma endregion

