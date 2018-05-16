/**************************************************************************
 Double List Template

 Author: Danilo Babin
 File name: "bdList.h"
***************************************************************************/



#ifndef BD_DOUBLE_LIST_DEF
	#define BD_DOUBLE_LIST_DEF


#include <typeinfo>
#include <cstdlib>
#include <sstream>

#include "bdArray.h"



template <class T> class bdListNode;
template <class T> class bdList;



//-------------------------------------------------------------------------------------------------------------------------------




template <class T>
class bdListNode
{
public:
	/// Element contained in the node.
	T m_element;
	/// Pointers to right and left node.
	bdListNode<T> *m_left, *m_right;
	/// Pointer to list that contains this node.
	bdList<T> *m_list;

	/// Constructor / Destructor
	bdListNode(bdList<T> *list);
	~bdListNode(){this->Reset();};
	
	/// Reset to the initial state (just after construction).
	void Reset();

	/// Get element contained in the node.
	T& GetElement(){return m_element;};

	/// Get Element to the right.
	bdListNode<T>* GetRight(){ return m_right; }

	/// Get Element to the left.
	bdListNode<T>* GetLeft(){ return m_left; }

	/// Get pointer to element contained in the node.
	T* GetElementPointer(){return (&m_element);};

	/// Get pointer to list to which the node belongs.
	bdList<T>* GetList(){return m_list;};

	/// Checks if it is the right/left end in the list it belongs to.
	int IsRightEnd();
	int IsLeftEnd();

	/// Print node type and header info.
	ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	ostream& PrintHeader(ostream &o);
};



template<class T>
void bdListNode<T>::Reset()
{
	//The list pointer remains the same - not reset!
	m_left = NULL; 
	m_right = NULL; 
	//m_element.~T();
}

template<class T>
bdListNode<T>::bdListNode(bdList<T> *list)
{	
	m_left = NULL; 
	m_right = NULL;
	m_list = list;
}

template<class T> 
int bdListNode<T>::IsRightEnd()
{
	if(this == m_list->GetRightEndNodePointer()) return 1;
	else return 0;
}
	
template<class T> 
int bdListNode<T>::IsLeftEnd()
{
	if(m_left == m_list->GetLeftEndNodePointer()) return 1;
	else return 0;
}

template<class T> 
ostream& bdListNode<T>::PrintHeader(ostream &o)
{
	PrintType(o)<<endl;
	o<<" - parent list address: "<<((void*)m_list)<<endl;
	o<<" - element address: "<<((void*)GetElementPointer())<<endl;
	o<<" - left node address: "<<((void*)m_left)<<endl;
	o<<" - right node address: "<<((void*)m_right)<<endl;
	return(o);
}



//-------------------------------------------------------------------------------------------------------------------------------

/// The bdListIterator should be used like this: 
/// for(iterator.SetLeftEnd(list); iterator.IsValid(); iterator.MoveToRight()) {...do something}


template<class T>
class bdListIterator
{
private:
	bdListNode<T> *m_node;

public:

	/// Constructor / Destructor.
	bdListIterator(){m_node = NULL;};
	~bdListIterator(){this->Reset();};

	/// Reset object to initial state (just after construction).
	void Reset(){m_node = NULL;};

	/// Set the iterator to given node, left or right of a list.
	void Set(bdListNode<T> *node){this->m_node = node;};
	void SetLeftEnd(bdList<T> &list){this->Set(list.GetLeftEndNodePointer());};
	void SetLeftEnd(bdList<T> *list){this->Set(list->GetLeftEndNodePointer());};
	void SetRightEnd(bdList<T> &list){this->Set(list.GetRightEndNodePointer());};
	void SetRightEnd(bdList<T> *list){this->Set(list->GetRightEndNodePointer());};

	/// Check if the iterator is valid.
	int IsValid();

	/// Move iterator in the set list to right/left with or without validity check.
	int MoveLeftWithValidityCheck();
	void MoveLeft(){m_node = m_node->m_left;};//Does not check validity of node pointer (assumes it is correct- use in FOR loop!)
	int MoveRightWithValidityCheck();
	void MoveRight(){m_node = m_node->m_right;};//Does not check validity of node pointer (assumes it is correct- use in FOR loop!)
	
	/// Get element or its pointer.
	int GetElementPointerWithValidityCheck(T **p_element);
	T* GetElementPointer(){return(m_node->GetElementPointer());};//Does not check validity of node pointer (assumes it is correct- use in FOR loop!)
	T& GetElement(){return(*(m_node->GetElementPointer()));};
	
	/// Get Node address.
	bdListNode<T>* GetNodeAddress(){return m_node;};
	bdListIterator<T>& operator =(bdListNode<T> *n);

	/// Print node type and header info.
	ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	ostream& PrintHeader(ostream &o);
};


template<class T> 
int bdListIterator<T>::IsValid()
{
	if(this->m_node==NULL) return 0;
	else return 1;
}

template<class T> 
int bdListIterator<T>::MoveLeftWithValidityCheck()
{
	if(this->m_node==NULL) return 0;
	else
	{
		m_node = m_node->left;
		return 1;
	}
}

template<class T> 
int bdListIterator<T>::MoveRightWithValidityCheck()
{
	if(this->m_node==NULL) return 0;
	else
	{
		m_node = m_node->right;
		return 1;
	}
}

template<class T> 
int bdListIterator<T>::GetElementPointerWithValidityCheck(T **p_element)
{
	if(this->m_node==NULL) 
	{
		*p_element = NULL;
		return 0;
	}
	else
	{
		*p_element = m_node->GetElementPointer();
		return 1;
	}
}


template<class T> 
bdListIterator<T>& bdListIterator<T>::operator =(bdListNode<T> *n)
{
	m_node = n;
}


template<class T> 
ostream& bdListIterator<T>::PrintHeader(ostream &o)
{
	PrintType(o)<<endl;
	o<<" - element address: "<<((void*)GetElementPointer())<<endl;
	o<<" - left node address: "<<((void*)m_node->m_left)<<endl;
	o<<" - right node address: "<<((void*)m_node->m_right)<<endl;
	return(o);
}





//-------------------------------------------------------------------------------------------------------------------------------





template <class T>
class bdList
{
private:
	/// Number of elements
	unsigned int m_noel;
	/// Pointer to left and right end nodes.
	bdListNode<T> *m_left_end, *m_right_end;
	
	/// Container node is the node that (might) contain this list (if there is a list of lists).
	bdListNode<T> *m_container_node;
	

public:

	/// Constructor / Destructor.
	bdList();
	bdList(bdListNode<T> *container_node);
	~bdList();

	/// Reset the object to intial state (just after constructor). The container node is NOT being reset!
	void Reset();

	/// Check if the list is empty. 
	int IsEmpty();

	/// Get number of elements in the list.
	unsigned int GetNumberOfElements() {return m_noel;};

	/// Get left end element, its pointer and node.
	T& GetLeftEnd();
	T* GetLeftEndPointerToElement();
	bdListNode<T>* GetLeftEndNodePointer(){return m_left_end;};

	/// Get right end element, its pointer and node.
	T& GetRightEnd();
	T* GetRightEndPointerToElement();
	bdListNode<T>* GetRightEndNodePointer(){return m_right_end;};

	/// Delete node, right and left end.
	int DeleteNode(bdListNode<T> *node);
	void DeleteLeftEnd();
	void DeleteRightEnd();

	/// Add element to list.
	void AddToLeftEnd(T &element);
	void AddToRightEnd(T &element);
	T* AddNewToRightEnd();
	T* AddNewToLeftEnd();
	int AddToLeftOfNode(bdListNode<T> *p_node, T &element);
	int AddToRightOfNode(bdListNode<T> *p_node, T &element);	

	/// Get Node pointer for given index or element. Return NULL if fail.
	bdListNode<T>* GetNodeWithIndex(unsigned int index);
	//bdListNode<T>* GetNodeWithElement(T &element);
	bdListNode<T>* GetNodeWithElementAddress(T *element_address);
    
    /// For this list creates a corresponding list with pointers to elements. If empty return fail 0, non-empty success 1.
    int CreateListOfPointersToElements(bdList<T*> &output_list_of_pointers);

	///// Efficient way to see if element is in ASCENDING SORTED list. Goes through the list and as soon as it finds that 
	///// current node element is higher than input element it returns fail (NULL). If element is found, returns node address.
	//bdListNode<T>* GetNodeWithElement_Ascending(T &element);

	///// Efficient way to get the first element greater than input (>) in ASCENDING SORTED list. Goes through the list and 
	///// as soon as it finds greater element, returns its pointer. If there is no greater element (end of list), return NULL.
	//bdListNode<T>* GetNodeWithElementGreaterThan_Ascending(T &element);

	///// Efficient way to get the first element greater than or equal to input (>=) in ASCENDING SORTED list. Goes through the list 
	///// and as soon as it finds greater element, returns its pointer. If there is no greater/equal element (end of list), return NULL.
	//bdListNode<T>* GetNodeWithElementGreaterOrEqual_Ascending(T &element);

	/// Add an input list to right end of this list.
	void AddListToRightEnd(bdList<T> &list);

	///// Check if the input elment is stored somewhere in the list. Element type must support operator '=='.
	//int HasElement(T &element);

	/// Get the index of input element or node. If element found return 1, else return 0.
	//int GetIndexOfElement(T &element, unsigned int &output_index);
	int GetIndexOfNode(bdListNode<T> *node, unsigned int &output_index);

	/// Puts the 'input list' to left end of this list, after the method input_list will be empty! 
	///  No copying is done, just pointer re-arrangement!
	void MergeInputListWithThisFromLeft(bdList<T> &input_list);

	/// Puts the 'input list' to right end of this list, after the method input_list will be empty! 
	///  No copying is done, just pointer re-arrangement!
	void MergeInputListWithThisFromRight(bdList<T> &input_list);

	/// Pointers are inverted: left end becomes right end, every left pointer becomes right pointer (and vice versa).
	///  DO NOT USE INSIDE LOOPS.
	void InvertOrder();

	/// Sort in ascending or descending order.
	//void SortAscending();
	//void SortDescending();
	int SortAscendingBySortingInputArray(bdArray<double> &input_array_to_sort);

	///Operators
	T& operator [](int r);
	bdList<T>& operator =(bdList<T> &l);

	/// Print info about type and header.
	ostream& PrintHeader(ostream &o);
	ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	friend std::ostream& operator <<(ostream &o, bdList<T> &l);
};


//template <class T>
//class bdListProcessing
//{
//public:
//
//	/// Constructor / Destructor.
//	bdListProcessing(){};
//	~bdListProcessing(){};
//
//	/// Efficient way to see if element is in ASCENDING SORTED list. Goes through the list and as soon as it finds that 
//	/// current node element is higher than input element it returns fail (NULL). If element is found, returns node address.
//	bdListNode<T>* GetNodeWithElement_Ascending(bdList<T> &list, T &element);
//
//	/// Efficient way to get the first element greater than input (>) in ASCENDING SORTED list. Goes through the list and 
//	/// as soon as it finds greater element, returns its pointer. If there is no greater element (end of list), return NULL.
//	bdListNode<T>* GetNodeWithElementGreaterThan_Ascending(bdList<T> &list, T &element);
//
//	/// Efficient way to get the first element greater than or equal to input (>=) in ASCENDING SORTED list. Goes through the list 
//	/// and as soon as it finds greater element, returns its pointer. If there is no greater/equal element (end of list), return NULL.
//	bdListNode<T>* GetNodeWithElementGreaterOrEqual_Ascending(bdList<T> &list, T &element);
//
//	/// Return pointer to node input with element. The element type must support '==' operator.
//	bdListNode<T>* GetNodeWithElement(bdList<T> &list, T &element);
//
//	/// Get index of element in the list. The element type must support '==' operator. If none found return 0, else 1.
//	int GetIndexOfElement(bdList<T> &list, T &element, unsigned int &output_index)
//
//	/// Check if the list has the input element. The element type must support '==' operator.
//	int HasElement(bdList<T> &list, T &element);
//
//	/// Sort elements in the list in the ascending order.
//	void SortAscending(bdList<T> &list);
//
//	/// Sort elements in the list in the descending order.
//	void SortDescending(bdList<T> &list);
//
//};


namespace bdListing
{
	/// Efficient way to see if element is in ASCENDING SORTED list. Goes through the list and as soon as it finds that 
	/// current node element is higher than input element it returns fail (NULL). If element is found, returns node address.
	template <typename T>
	bdListNode<T>* GetNodeWithElement_Ascending(bdList<T> &list, T &element);

	/// Efficient way to get the first element greater than input (>) in ASCENDING SORTED list. Goes through the list and 
	/// as soon as it finds greater element, returns its pointer. If there is no greater element (end of list), return NULL.
	template <typename T>
	bdListNode<T>* GetNodeWithElementGreaterThan_Ascending(bdList<T> &list, T &element);

	/// Efficient way to get the first element greater than or equal to input (>=) in ASCENDING SORTED list. Goes through the list 
	/// and as soon as it finds greater element, returns its pointer. If there is no greater/equal element (end of list), return NULL.
	template <typename T>
	bdListNode<T>* GetNodeWithElementGreaterOrEqual_Ascending(bdList<T> &list, T &element);

	/// Return pointer to node input with element. The element type must support '==' operator.
	template <typename T>
	bdListNode<T>* GetNodeWithElement(bdList<T> &list, T &element);

	/// Get index of element in the list. The element type must support '==' operator. If none found return 0, else 1.
	template <typename T>
	int GetIndexOfElement(bdList<T> &list, T &element, unsigned int &output_index);

	/// Check if the list has the input element. The element type must support '==' operator.
	template <typename T>
	int HasElement(bdList<T> &list, T &element);

	/// Sort elements in the list in the ascending order.
	template <typename T>
	void SortAscending(bdList<T> &list);

	/// Sort elements in the list in the descending order.
	template <typename T>
	void SortDescending(bdList<T> &list);

};



template<class T> class bdList;  // pre-declare the template class itself
//template<class T> std::ostream& operator<< (ostream &o, DoubleList<T> &l); 



template <class T>
bdList<T>::bdList()
{
    m_noel = 0;
	m_left_end = m_right_end = NULL;
	m_container_node = NULL;
}


template <class T>
bdList<T>::bdList(bdListNode<T> *container_node)
{
    m_noel = 0;
	m_left_end = m_right_end = NULL;
	m_container_node = container_node;
}


template <class T>
bdList<T>::~bdList()
{
	this->Reset();
}


template<class T>
void bdList<T>::Reset()
{
	while(!this->IsEmpty()) this->DeleteRightEnd();
	m_left_end = m_right_end = NULL;
	m_noel = 0;
	// The container node is NOT being reset!
}


template <class T>
int bdList<T>::IsEmpty()
{
    if(m_noel<=0) return 1;
	else return 0;
}


template <class T>
T& bdList<T>::GetLeftEnd()
{
	return (m_left_end->GetElement());
}


template <class T>
T& bdList<T>::GetRightEnd()
{
	return (m_right_end->GetElement());
}


template <class T>
T* bdList<T>::GetLeftEndPointerToElement()
{
	return (m_left_end->GetElementPointer());
}


template <class T>
T* bdList<T>::GetRightEndPointerToElement()
{
	return (m_right_end->GetElementPointer());
}


template<class T>
int bdList<T>::DeleteNode(bdListNode<T> *node)
{
	if(node==NULL) return 0;

	//If there is only 1 element in the list
	if(m_noel == 1)
	{
		delete m_right_end;
		m_left_end = m_right_end = NULL;
		m_noel = 0;
		return 1;
	}

	//If the current_node is the right_end
	if(node==m_right_end)
	{
		this->DeleteRightEnd();
		return 1;
	}

	//If the current_node is the left_end
	if(node==m_left_end)
	{
		this->DeleteLeftEnd();
		return 1;
	}

	node->m_right->m_left = node->m_left;
	node->m_left->m_right = node->m_right;
	delete node;
	m_noel--;

	return 1;
}


template <class T>
void bdList<T>::DeleteLeftEnd()
{
	if(this->IsEmpty()) return;

	if(m_noel == 1)
	{		
	    delete m_left_end;
	    m_left_end = m_right_end = NULL;
	    m_noel = 0;
	}
	else
	{
		bdListNode<T> *del;
	    del = m_left_end;
	    m_left_end = m_left_end->m_right;
	    m_left_end->m_left = NULL;		
	    delete del;
	    m_noel--;
	}
}


template <class T>
void bdList<T>::DeleteRightEnd()
{
	if(this->IsEmpty()) return;

	if(m_noel == 1)
	{
		delete m_right_end;
	    m_left_end = m_right_end = NULL;
	    m_noel = 0;
	}
	else
	{		
		bdListNode<T> *del;
	    del = m_right_end;
	    m_right_end = m_right_end->m_left;
	    m_right_end->m_right = NULL;
	    delete del;
	    m_noel--;
	}
}


template<class T>
void bdList<T>::AddToRightEnd(T &element)
{
	T *nel = this->AddNewToRightEnd();
	(*nel) = element;
}


template<class T>
T* bdList<T>::AddNewToRightEnd()
{
    bdListNode<T> *nel = new bdListNode<T>(this);
	if(this->IsEmpty())
    {
		nel->m_left = NULL;
		nel->m_right = NULL;
		m_left_end = m_right_end = nel;
		m_noel++;
    }
    else
    {
		nel->m_right = NULL;
		nel->m_left = m_right_end;
		m_right_end->m_right = nel;
		m_right_end = nel;
		m_noel++;
    }
	return (nel->GetElementPointer());
}


template<class T>
void bdList<T>::AddToLeftEnd(T &element)
{
	T *nel = this->AddNewToLeftEnd();
	(*nel) = element;
}


template<class T>
T* bdList<T>::AddNewToLeftEnd()
{
    bdListNode<T> *nel = new bdListNode<T>(this);
	if(this->IsEmpty())
    {
		nel->m_left = NULL;
		nel->m_right = NULL;
		m_left_end = m_right_end = nel;
		m_noel++;
    }
    else
    {
		nel->m_left = NULL;
		nel->m_right = m_left_end;
		m_left_end->m_left = nel;
		m_left_end = nel;
		m_noel++;
    }
	return (nel->GetElementPointer());
}


template<class T>
int bdList<T>::AddToLeftOfNode(bdListNode<T> *node, T &element)
{
	if(!node) return 0;
	if(this->IsEmpty()) return 0;

	if(node==m_left_end){ this->AddToLeftEnd(element); }
	else 
	{
		bdListNode<T> *nel = new bdListNode<T>(this);
		nel->GetElement() = element;

		nel->m_left = node->m_left;
		nel->m_right = node;
		node->m_left->m_right = nel;
		node->m_left = nel;
		m_noel++;
	}
	return 1;
}


template<class T>
int bdList<T>::AddToRightOfNode(bdListNode<T> *node, T &element)
{
	if(!node) return 0;
	if(this->IsEmpty()) return 0;

	if(node==m_right_end){ this->AddToRightEnd(element); }
	else 
	{
		bdListNode<T> *nel = new bdListNode<T>(this);
		nel->GetElement() = element;
		nel->m_right = node->m_right;
		nel->m_left = node;
		node->m_right->m_left = nel;
		node->m_right = nel;
		m_noel++;
	}
	return 1;
}


template<class T>
bdListNode<T>* bdList<T>::GetNodeWithIndex(unsigned int index)
{
	if(index>=m_noel) return NULL;
	bdListNode<T> *node = m_left_end;
	for(unsigned int i=0; i<index; i++) node = node->m_right;
	return node;
}


//template<class T>
//bdListNode<T>* bdList<T>::GetNodeWithElement(T &element)
//{
//	if(this->IsEmpty()) return NULL;
//	bdListIterator<T> it;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( it.GetElement() == element ) 
//		{
//			return it.GetNodeAddress();
//		}
//	}
//	return NULL;
//}


template<class T>
bdListNode<T>* bdList<T>::GetNodeWithElementAddress(T *element_address)
{
	if(this->IsEmpty()) return NULL;
	bdListIterator<T> it;
	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
	{
		if( it.GetElementPointer() == element_address ) 
		{
			return it.GetNodeAddress();
		}
	}
	return NULL;
}


template<class T>
int bdList<T>::CreateListOfPointersToElements(bdList<T*> &output_list_of_pointers)
{
    if(this->IsEmpty()) return 0;
    output_list_of_pointers.Reset();
    
    bdListIterator<T> it;
    for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
    {
        T *p = it.GetElementPointer();
        output_list_of_pointers.AddToRightEnd(p);
    }
    
    return 1;
}

//template<class T>
//bdListNode<T>* bdList<T>::GetNodeWithElement_Ascending(T &element)
//{
//	if(this->IsEmpty()) return NULL;
//	bdListIterator<T> it;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( it.GetElement() == element ) 
//		{
//			return it.GetNodeAddress();
//		}
//		else
//		{
//			if(it.GetElement() > element) return NULL;
//		}
//	}
//	return NULL;
//}
//
//
//template<class T>
//bdListNode<T>* bdList<T>::GetNodeWithElementGreaterThan_Ascending(T &element)
//{
//	if(this->IsEmpty()) return NULL;
//	bdListIterator<T> it;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( it.GetElement() > element ) { return it.GetNodeAddress(); }
//	}
//	return NULL;
//}
//
//
//template<class T>
//bdListNode<T>* bdList<T>::GetNodeWithElementGreaterOrEqual_Ascending(T &element)
//{
//	if(this->IsEmpty()) return NULL;
//	bdListIterator<T> it;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( it.GetElement() >= element ) { return it.GetNodeAddress(); }
//	}
//	return NULL;
//}


//template<class T>
//int bdList<T>::HasElement(T &element)
//{
//	if(this->IsEmpty()) return 0;
//	bdListIterator<T> it;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( (*(it.GetElementPointer()))==element ) return 1;
//	}
//	
//	return 0;
//}


//template<class T>
//int bdList<T>::GetIndexOfElement(T &element, unsigned int &output_index)
//{
//	if(this->IsEmpty()) return 0;
//	bdListIterator<T> it;
//	output_index = 0;
//	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
//	{
//		if( (*(it.GetElementPointer()))==element ) return 1;
//		output_index++;
//	}
//	return 0;
//}


template<class T>
int bdList<T>::GetIndexOfNode(bdListNode<T> *node, unsigned int &output_index)
{
	if(this->IsEmpty()) return 0;
	bdListIterator<T> it;
	output_index = 0;
	for(it.SetLeftEnd(this); it.IsValid(); it.MoveRight())
	{
		if( it.GetNodeAddress()==node ) return 1;
		output_index++;
	}
	return 0;
}


template<class T>
void bdList<T>::AddListToRightEnd(bdList<T> &list)
{
	if(list.IsEmpty()) return;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		this->AddToRightEnd(it.GetElement());
	}
}


template<class T>
void bdList<T>::MergeInputListWithThisFromLeft(bdList<T> &input_list)
{
	//If the input list is empty do nothing
	if(input_list.IsEmpty()) return;

	// Set the pointer of nodes in the input list to point to this list 
	bdListIterator<T> it;
	for(it.SetLeftEnd(input_list); it.IsValid(); it.MoveRight())
	{ it.GetNodeAddress()->m_list = this; }

	//If this list is empty
	if(this->IsEmpty())
	{
		//Set the pointers of this list to point to the data of input_list
		m_right_end = input_list.m_right_end;
		m_left_end = input_list.m_left_end;
		m_noel = input_list.m_noel;
	}
	//If this list is not empty
	else
	{
		//Set the pointers of this list to point to the data of input_list
		m_left_end->m_left = input_list.m_right_end;
		input_list.m_right_end->m_right = m_left_end;
		m_left_end = input_list.m_left_end;
		m_noel += input_list.m_noel;//Set number of elements (it is sum of number of elements of both lists)
	}

	//Detach the input list from its data
	input_list.m_left_end = NULL;
	input_list.m_right_end = NULL;
	input_list.m_noel = 0;
}


template<class T>
void bdList<T>::MergeInputListWithThisFromRight(bdList<T> &input_list)
{
	// If the input list is empty do nothing
	if(input_list.IsEmpty()) return;

	// Set the pointer of nodes in the input list to point to this list 
	bdListIterator<T> it;
	for(it.SetLeftEnd(input_list); it.IsValid(); it.MoveRight())
	{ it.GetNodeAddress()->m_list = this; }

	// If this list is empty
	if(this->IsEmpty())
	{
		//Set the pointers of this list to point to the data of input_list
		m_right_end = input_list.m_right_end;
		m_left_end = input_list.m_left_end;
		m_noel = input_list.m_noel;
	}
	// If this list is not empty
	else
	{
		//Set the pointers of this list to point to the data of input_list
		m_right_end->m_right = input_list.m_left_end;
		input_list.m_left_end->m_left = m_right_end;
		m_right_end = input_list.m_right_end;
		m_noel += input_list.m_noel;//Set number of elements (it is sum of number of elements of both lists)
	}

	// Detach the input list from its data
	input_list.m_left_end = NULL;
	input_list.m_right_end = NULL;
	input_list.m_noel = 0;
}


template<class T>
void bdList<T>::InvertOrder()
{
	if(this->IsEmpty()) return;

	bdListNode<T> *temp, *temp2, *current;

	// Invert all left and right pointers in nodes
	current = this->GetLeftEndNodePointer();
	while(current!=m_right_end)
	{
		temp = current;
		current = current->m_right;
		temp2 = temp->m_left;
		temp->m_left = temp->m_right;
		temp->m_right = temp2;
	}
	// At this point we have to do it also for the last node
	temp2 = current->m_left;
	current->m_left = current->m_right;
	current->m_right = temp2;

	// Invert left end and right end
	temp = m_left_end;
	m_left_end = m_right_end;
	m_right_end = temp;
}


//template<class T>
//void bdList<T>::SortAscending()
//{
//	if(m_noel<=1) return;
//
//	// Indicator if any changes were made in the while loop
//	int is_change_made = 1;
//
//	bdListNode<T> *node1, *node2;
//	T temp;
//
//	// Loop while there are changes in order
//	while(is_change_made)
//	{
//		// Reset indicator, initialize pointers
//		is_change_made = 0;
//		node1 = m_left_end;
//		node2 = m_left_end->m_right;
//
//		for(int i=0; i<m_noel-1; i++)// can be also: while(node2!=NULL), but 'for' loop with counter is safer here 
//		{
//			if( (node1->GetElement()) > (node2->GetElement()) )
//			{
//				is_change_made = 1;
//				temp = node2->GetElement();
//				node2->GetElement() = node1->GetElement();
//				node1->GetElement() = temp;
//			}
//			node1 = node1->m_right;
//			node2 = node2->m_right;
//		}
//	}
//}
//
//
//template<class T>
//void bdList<T>::SortDescending()
//{
//	if(m_noel<=1) return;
//
//	// Indicator if any changes were made in the while loop
//	int is_change_made = 1;
//
//	bdListNode<T> *node1, *node2;
//	T temp;
//
//	// Loop while there are changes in order
//	while(is_change_made)
//	{
//		// Reset indicator, initialize pointers
//		is_change_made = 0;
//		node1 = m_left_end;
//		node2 = m_left_end->m_right;
//
//		for(int i=0; i<m_noel-1; i++)// can be also: while(node2!=NULL), but 'for' loop with counter is safer here 
//		{
//			if( (node1->GetElement()) < (node2->GetElement()) )
//			{
//				is_change_made = 1;
//				temp = node2->GetElement();
//				node2->GetElement() = node1->GetElement();
//				node1->GetElement() = temp;
//			}
//			node1 = node1->m_right;
//			node2 = node2->m_right;
//		}
//	}
//}


template <class T>
int bdList<T>::SortAscendingBySortingInputArray(bdArray<double> &input_array_to_sort)
{
	if(m_noel<=0) return 0;
	if(input_array_to_sort.GetNumberOfElements()!=m_noel) return 0;
	if(m_noel==1) return 1;

	// Indicator if any changes were made in the while loop
	int is_change_made = 1;

	bdListNode<T> *p_node1, *p_node2, node_temp(NULL), *p_node_iterator;

	// Loop while there are changes in order
	while(is_change_made)
	{
		// Reset indicator, initialize pointers
		is_change_made = 0;
		p_node_iterator = m_left_end;

		for(unsigned int i=0; i<m_noel-1; i++)// can be also: while(node2!=NULL), but 'for' loop with counter is safer here 
		{
			if(input_array_to_sort[i] > input_array_to_sort[i+1])
			{
				is_change_made = 1;

				//Make the change in the array
				double temp = input_array_to_sort[i];
				input_array_to_sort[i] = input_array_to_sort[i+1];
				input_array_to_sort[i+1] = temp;

				//Make the change in the list
				p_node1 = p_node_iterator;
				p_node2 = p_node_iterator->m_right;

				node_temp.m_left = p_node1->m_left;
				node_temp.m_right = p_node1->m_right;

				p_node1->m_left = p_node2;
				p_node1->m_right = p_node2->m_right;
				
				if(node_temp.m_left!=NULL) node_temp.m_left->m_right = p_node2;
				else this->m_left_end = p_node2;
				if(p_node1->m_right!=NULL) p_node1->m_right->m_left = p_node1;
				else this->m_right_end = p_node1;
				p_node2->m_left = node_temp.m_left;
				p_node2->m_right = p_node1;

				//Set the iterator
				p_node_iterator = p_node2;//Set the iterator pointer (p_node2 is now left from p_node1, so that's why we use p_node2).
			}
			p_node_iterator = p_node_iterator->m_right;
		}
	}
	return 1;
}


template <class T>
T& bdList<T>::operator [](int r)
{ 
	bdListNode<T> *node;
	node = m_left_end;
	for(int i=0; i<r; i++) node = node->m_right;
	return (node->GetElement());
}

template <class T>
bdList<T>& bdList<T>::operator =(bdList<T> &l)
{
    if (&l==this) return *this;
	this->Reset();
	bdListIterator<T> it;
	for(it.SetLeftEnd(l); it.IsValid(); it.MoveRight()) { this->AddToRightEnd(it.GetElement()); }
    return *this;
}

template <class T>
ostream& bdList<T>::PrintHeader(ostream &o)
{
	PrintType(o)<<endl;
	o<<" - parent container address: "<<((void*)m_container_node)<<endl;
	o<<" - number of elements: "<<this->GetNumberOfElements()<<endl;
	o<<" - left end address: "<<((void*)m_left_end)<<endl;
	o<<" - right end address: "<<((void*)m_right_end)<<endl;
	return(o);
}

template <class T>
std::ostream& operator <<(ostream &o, bdList<T> &l)
{
	bdListIterator<T> it;
	for(it.SetLeftEnd(l); it.IsValid(); it.MoveRight())
	{
		o<<it.GetElement()<<" ";
	}
	o<<endl;
    return o;
}




//---------------------------------------------------------------------------------------------------------------------




template<class T>
bdListNode<T>* bdListing::GetNodeWithElement_Ascending(bdList<T> &list, T &element)
{
	if(list.IsEmpty()) return NULL;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( it.GetElement() == element ) 
		{
			return it.GetNodeAddress();
		}
		else
		{
			if(it.GetElement() > element) return NULL;
		}
	}
	return NULL;
}


template<class T>
bdListNode<T>* bdListing::GetNodeWithElementGreaterThan_Ascending(bdList<T> &list, T &element)
{
	if(list.IsEmpty()) return NULL;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( it.GetElement() > element ) { return it.GetNodeAddress(); }
	}
	return NULL;
}


template<class T>
bdListNode<T>* bdListing::GetNodeWithElementGreaterOrEqual_Ascending(bdList<T> &list, T &element)
{
	if(list.IsEmpty()) return NULL;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( it.GetElement() >= element ) { return it.GetNodeAddress(); }
	}
	return NULL;
}


template<class T>
bdListNode<T>* bdListing::GetNodeWithElement(bdList<T> &list, T &element)
{
	if(list.IsEmpty()) return NULL;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( it.GetElement() == element ) 
		{
			return it.GetNodeAddress();
		}
	}
	return NULL;
}


template<class T>
int bdListing::GetIndexOfElement(bdList<T> &list, T &element, unsigned int &output_index)
{
	if(list.IsEmpty()) return 0;
	bdListIterator<T> it;
	output_index = 0;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( it.GetElement() == element ) return 1;
		output_index++;
	}
	return 0;
}



template<class T>
int bdListing::HasElement(bdList<T> &list, T &element)
{
	if(list.IsEmpty()) return 0;
	bdListIterator<T> it;
	for(it.SetLeftEnd(list); it.IsValid(); it.MoveRight())
	{
		if( (*(it.GetElementPointer()))==element ) return 1;
	}
	return 0;
}


template<class T>
void bdListing::SortAscending(bdList<T> &list)
{
	if(list.GetNumberOfElements()<=1) return;

	// Indicator if any changes were made in the while loop
	int is_change_made = 1;

	bdListNode<T> *node1, *node2;
	T temp;

	// Loop while there are changes in order
	while(is_change_made)
	{
		// Reset indicator, initialize pointers
		is_change_made = 0;
		node1 = list.GetLeftEndNodePointer();
		node2 = list.GetLeftEndNodePointer()->GetRight();

		for(unsigned int i=0; i<list.GetNumberOfElements()-1; i++)// can be also: while(node2!=NULL), but 'for' loop with counter is safer here 
		{
			if( (node1->GetElement()) > (node2->GetElement()) )
			{
				is_change_made = 1;
				temp = node2->GetElement();
				node2->GetElement() = node1->GetElement();
				node1->GetElement() = temp;
			}
			node1 = node1->GetRight();
			node2 = node2->GetRight();
		}
	}
}


template<class T>
void bdListing::SortDescending(bdList<T> &list)
{
	if(list.GetNumberOfElements()<=1) return;

	// Indicator if any changes were made in the while loop
	int is_change_made = 1;

	bdListNode<T> *node1, *node2;
	T temp;

	// Loop while there are changes in order
	while(is_change_made)
	{
		// Reset indicator, initialize pointers
		is_change_made = 0;
		node1 = list.GetLeftEndNodePointer();
		node2 = list.GetLeftEndNodePointer()->GetRight();

		for(unsigned int i=0; i<list.GetNumberOfElements()-1; i++)// can be also: while(node2!=NULL), but 'for' loop with counter is safer here 
		{
			if( (node1->GetElement()) < (node2->GetElement()) )
			{
				is_change_made = 1;
				temp = node2->GetElement();
				node2->GetElement() = node1->GetElement();
				node1->GetElement() = temp;
			}
			node1 = node1->GetRight();
			node2 = node2->GetRight();
		}
	}
}


#endif