#include "Event.h"
//copied from my previous version OI 15.3

#pragma region ctors
Event::Event(bool indel_flag, int ind, int len)
{
	this->index = ind;
	this->length = len;
	this->indel_flag = indel_flag;
}
Event::Event()
{
}
#pragma endregion
#pragma region dtrs
Event::~Event()
{

}
#pragma endregion
#pragma region funcs
bool Event::getIndel_Flag() {
	return this->indel_flag;
}
int Event::getIndex() {
	return this->index;
}
int Event::getLength() {
	return this->length;
}
void Event::setIndel_Flag(bool indel) {
	this->indel_flag = indel;
}
void Event::setIndex(int ind) {
	this->index = ind;
	return;
}
void Event::setLength(int len) {
	this->length = len;
	return;
}
#pragma endregion



