#pragma once
class Event
{
	//copied from my previous version OI 15.3

private:
	bool indel_flag;
	int index;
	int length;
public:
	Event();
	Event(bool indel_flag, int ind, int len);

	~Event();
	
	void setIndex(int ind);
	void setLength(int len);
	void setIndel_Flag(bool indel);

	int getLength();
	int getIndex();
	bool getIndel_Flag();

};

