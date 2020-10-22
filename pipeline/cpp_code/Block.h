#pragma once
//copied from my previous version OI 15.3

class Block
{
private:
	int index;
	int length;
public:
	Block();
	Block(int index, int length);
	~Block();
	bool operator ==(Block block);
	void setIndex(int index);
	void setLength(int length);
	int getLength();
	int getIndex();

};

