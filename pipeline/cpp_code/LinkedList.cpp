#include "LinkedList.h"

using namespace std;
#pragma region ctors
LinkedList::LinkedList()
{
	this->head = NULL;
	this->tail = NULL;
	this->length = 0;
}
#pragma endregion
#pragma region dtors
LinkedList::~LinkedList() {
}
void LinkedList::delete_LL()
{
    Node* next = head;

   while (next!=NULL) {              // iterate over all elements
       Node* deleteMe = next;
       next = (*next).getNext();     // save pointer to the next element
       delete deleteMe;       // delete the current entry
   }
}

#pragma endregion
#pragma region FUNCS
void LinkedList::print()
{
    Node* pos = head;
    while (pos != tail) {
        cout << (*pos).getData() << " ";
        pos = (*pos).getNext();
    }

    cout <<(*tail).getData()<< endl;
}
size_t LinkedList::getLength()
{
    return length;
}
Node* LinkedList::getHead()
{
    return head;
}
Node* LinkedList::getTail()
{
    return tail;
}

void LinkedList::add_list_to_list(LinkedList& ll, Node* insertafter)
{
    //we assume ll isnt empty  and that insert after is in the node//18.3
    if (head==NULL)
    {//extreme condition- if the linkedlist is empty //OI 18.3
        head = ll.head;
        tail = ll.tail;
    }
    else
    {
        if (insertafter == NULL)// it means that we need to insert the list at the begining... \\OI 15.3
        {
            Node* temp = head;
            head = ll.head;
            ll.tail->setNext(temp);

        }
        else// now we insert after the node "insertafter" \\OI 15.3
        {
            Node* temp = insertafter->getNext();
            insertafter->setNext(ll.head);
            ll.tail->setNext(temp);
            if ((*insertafter).getData() == (*tail).getData())// small change- we now compare the values(data)
            {//NEW addition- updating the tail if needed //16.3 OI
                tail = ll.tail;
            }
        }
      
    }
    length += ll.length;
}
void LinkedList::addnode_to_list(Node* insertnode, Node* insertafter)
{
    //inserts **after the "insertafternode"
    if (head == NULL)// change 18.3 //OI 
    {// extreme cond- if the linked list is empty
        head = insertnode;
        tail = insertnode;
    }
    else
    {
        if (insertafter == NULL)// it means that we need to insert the list at the begining... \\OI 15.3
        {

            Node* temp = head;
            head = insertnode;
            insertnode->setNext(temp);//adding to the length
        }
        else// now we insert after the node "insertafter" \\OI 15.3
        {
            Node* temp = (*insertafter).getNext();
            (*insertafter).setNext(insertnode);
            if ((*insertafter).getData() == (*tail).getData())
            {//NEW addition- updating the tail if needed //16.3 OI
                (*insertnode).setNext(NULL);
                tail = insertnode;
                
            }
            else
            {//change of syntax and order of the function.. //OI 18.3
               (*insertnode).setNext(temp);
            }
        }
    }
    length += 1;//adding to the length
}
#pragma endregion


