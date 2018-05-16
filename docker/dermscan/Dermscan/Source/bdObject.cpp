/**************************************************************************
 Object class as parent of all image analysis classes later on.
 
 Progress Counter class used in Function Object and Object classes to
 track the progress of method execution.

 Author: Danilo Babin
 File name: "bdObject.h"
***************************************************************************/



#include "bdObject.h"



bdProgressCounter::bdProgressCounter()
{
    m_value = 0;
    m_minimum = 0;
    m_maximum = 100;
    m_object = NULL;
}


bdProgressCounter::~bdProgressCounter()
{
    m_value = 0;
    m_minimum = 0;
    m_maximum = 100;
    if(m_object) this->DisconnectObject();
    m_object = NULL;
}


//void bdProgressCounter::SetObject(bdObject *object)
//{
//	//If the progress counter was being used by another object, disconnect them.
//	if(m_object)
//	{
//		if(object!=m_object) m_object->m_progress_counter = NULL;
//	}
//	m_object = object;
//}



void bdProgressCounter::DisconnectObject()
{
    //If the progress counter was being used by another object, disconnect them.
    if(m_object)
    {
        m_object->m_progress_counter = NULL;
        m_object = NULL;
    }
}



void bdProgressCounter::Set(int minimum, int maximum, int value)
{
    m_value = value;
    m_minimum = minimum;
    m_maximum = maximum;
}


void bdProgressCounter::SetMinimum(int min)
{
    m_minimum = min;
}


void bdProgressCounter::SetMaximum(int max)
{
    m_maximum = max;
}


void bdProgressCounter::SetValue(int value)
{
    m_value = value;
}


void bdProgressCounter::SetValueRelative(int relative_value)
{
    if(relative_value<0 || relative_value>100) m_value = m_minimum;
    m_value = m_minimum + (relative_value*(m_maximum-m_minimum))/100;
}


int bdProgressCounter::GetMinimum()
{
    return m_minimum;
}


int bdProgressCounter::GetMaximum()
{
    return m_maximum;
}


int bdProgressCounter::GetValue()
{
    return m_value;
}





//------------------------------------------------------------------------------------------------------





bdObject::bdObject()
{
	this->SetupObject();
}


bdObject::~bdObject()
{
	this->Reset();
}


void bdObject::SetupObject()
{
	m_progress_counter = NULL;
}


void bdObject::Reset()
{
	this->DisconnectProgressCounter();
}


bdString& bdObject::GetTypeName()
{
    return m_type_name;
}


int bdObject::ConnectProgressCounter(bdProgressCounter *progress_counter)
{
	if(!progress_counter) return 0;

	//If the progress counter is already in use...
	if(progress_counter->m_object != NULL)
	{
		//... by another Object, deny access.
		if(progress_counter->m_object != this) return 0;
		//... other wise, just signal that it is already conneted to this object.
		else return 1;
	}
	//If the progress counter is not in use...
	else
	{
		//... connect it to this object.
		progress_counter->m_object = this;
		m_progress_counter = progress_counter;
		return 1;
	}
}


int bdObject::DisconnectProgressCounter()
{
	if(!m_progress_counter) return 1;

	//If the progress counter is in use...
	if(m_progress_counter->m_object != NULL)
	{
		//... by another Object, deny access.
		if(m_progress_counter->m_object != this) return 0;
		//... other wise disconnect from this object
		else 	
		{
			m_progress_counter->m_object = NULL;
			m_progress_counter = NULL;
			return 1;
		}
	}
	//If the progress counter is not in use return success.
	else return 1;
}



//void bdObject::SetProgressCounter(bdProgressCounter *progress_counter)
//{
//	m_progress_counter = progress_counter;
//}


void bdObject::SetProgressCounterProperties(int min, int max, int value)
{
	if(m_progress_counter) m_progress_counter->Set(min,max,value);
//	m_progress_counter.Set(min,max,value);
}


void bdObject::SetProgressCounterRelativeValue(int relative_value)
{
	if(m_progress_counter) m_progress_counter->SetValueRelative(relative_value);
//	m_progress_counter.SetValueRelative(relative_value);
}


//void bdObject::RemoveProgressCounter()
//{
//	m_progress_counter = NULL;
//}


//-----------------------------------------------------------------------------------------------------------------------




bdString bdDataObjectTag::GetValueForGivenIdentifier(bdString identifier)
{
    bdString value; value.Assign("");
  
    //cout<<"identifier="<<identifier.C_String();
    
    
//    for(int i=0; i<m_identifiers.GetNumberOfElements(); i++)
//    {
//        cout<<"|"<<m_identifiers[i].C_String()<<"|";
//    }
    
    unsigned int index;
    if(bdListing::GetIndexOfElement(m_identifiers, identifier, index))
    {
        value.Assign(m_values[index]);
    }
    
    return value;
}


bdString bdDataObjectTag::GetValueForGivenIdentifier(const char *identifier)
{
    bdString id; id.Assign(identifier);
    return this->GetValueForGivenIdentifier(id);
}


int bdDataObjectTag::GetValueAsIntegerForGivenIdentifier(bdString identifier, int &output_integer)
{
    bdString bds = GetValueForGivenIdentifier(identifier);
    if(bds.IsEmpty()) return 0;
    if(bds.StringToInt(output_integer)) return 1;
    return 0;
}


int bdDataObjectTag::GetValueAsIntegerForGivenIdentifier(const char *identifier, int &output_integer)
{
    bdString id; id.Assign(identifier);
    return this->GetValueAsIntegerForGivenIdentifier(id, output_integer);
}


int bdDataObjectTag::GetValueAsDoubleForGivenIdentifier(bdString identifier, double &output_double)
{
    bdString bds = GetValueForGivenIdentifier(identifier);
    if(bds.IsEmpty()) return 0;
    if(bds.StringToDouble(output_double)) return 1;
    return 0;
}


int bdDataObjectTag::GetValueAsDoubleForGivenIdentifier(const char *identifier, double &output_double)
{
    bdString id; id.Assign(identifier);
    return this->GetValueAsDoubleForGivenIdentifier(id, output_double);
}


void bdDataObjectTag::AddIdentifierAndValue(bdString identifier, bdString value)
{
    m_identifiers.AddToRightEnd(identifier);
    m_values.AddToRightEnd(value);
}


void bdDataObjectTag::AddIdentifierAndValue(const char *identifier, const char *value)
{
    bdString id; id.Assign(identifier);
    bdString v; v.Assign(value);
    this->AddIdentifierAndValue(id, v);
}


void bdDataObjectTag::AddDataObjectReference(bdString data_object_name, bdString value)
{
    bdString identifier; identifier.Assign("*"); identifier.Append(data_object_name);
    this->AddIdentifierAndValue(identifier, value);
}


void bdDataObjectTag::AddDataObjectReference(const char *data_object_name, const char *value)
{
    bdString identifier; identifier.Assign("*"); identifier.Append(data_object_name);
    bdString v; v.Assign(value);
    this->AddIdentifierAndValue(identifier, v);
}


int bdDataObjectTag::GetValueForGivenDataObjectName(bdString data_object_name, bdString &output_value)
{
    bdString identifier; identifier.Assign("*"); identifier.Append(data_object_name);
    unsigned int index;
    if(bdListing::GetIndexOfElement(m_identifiers, identifier, index))
    {
        output_value.Assign(m_values[index]);
        return 1;
    }
    return 0;
}


int bdDataObjectTag::GetValueForGivenDataObjectName(const char *data_object_name, bdString &output_value)
{
    bdString id; id.Assign(data_object_name);
    return (this->GetValueForGivenDataObjectName(id, output_value));
}


int bdDataObjectTag::HasDataObjectReference(bdString data_object_name, bdString value)
{
    bdString identifier; identifier.Assign("*"); identifier.Append(data_object_name);
    
    bdListIterator<bdString> iti, itv;
    for(iti.SetLeftEnd(m_identifiers), itv.SetLeftEnd(m_values); iti.IsValid() && itv.IsValid(); iti.MoveRight(), itv.MoveRight())
    {
        if(iti.GetElement()==identifier && itv.GetElement()==value) return 1;
    }
    return 0;
}


int bdDataObjectTag::HasDataObjectReference(const char *data_object_name, const char *value)
{
    bdString obj_name, v; obj_name.Assign(data_object_name); v.Assign(value);
    return (this->HasDataObjectReference(obj_name, v));
}


int bdDataObjectTag::HasDataObjectReference(bdString data_object_name)
{
    bdString identifier; identifier.Assign("*"); identifier.Append(data_object_name);
    
    bdListIterator<bdString> iti;
    for(iti.SetLeftEnd(m_identifiers); iti.IsValid(); iti.MoveRight())
    {
        if(iti.GetElement()==identifier) return 1;
    }
    return 0;
}


int bdDataObjectTag::HasDataObjectReference(const char *data_object_name)
{
    bdString obj_name; obj_name.Assign(data_object_name);
    return (this->HasDataObjectReference(obj_name));
}


bdDataObjectTag& bdDataObjectTag::operator =(bdDataObjectTag &tag)
{
    if(&tag==this) return *this;
    m_namespace = tag.m_namespace;
    m_class = tag.m_class;
    m_object = tag.m_object;
    m_identifiers  = tag.m_identifiers;
    m_values = tag.m_values;
    return *this;
}




//-----------------------------------------------------------------------------------------------------------------------


bdDataObject::bdDataObject()
{
	this->SetupObject();
	m_is_locked = 0;
    this->m_color_RGBA[3] = 0;//set to no color.
    this->m_has_exclusive_visibility = 0;// set to view multiple object simultaniously.
    this->GetObjectName().Assign(" ");
    this->AddTag("bd", "bdDataObject", " ");// basic tag that is always present when the object is created.
}


bdDataObject::~bdDataObject()
{
	this->Reset();
	m_is_locked = 0;
	m_history.Clear();
}


bdString& bdDataObject::GetHistory()
{
	return m_history;
}


bdString& bdDataObject::GetObjectName()
{
    return m_object_name;
}


int bdDataObject::HasExclusiveVisibility()
{
    return m_has_exclusive_visibility;
}


void bdDataObject::SetExclusiveVisibility(int is_visibility_exclusive)
{
    m_has_exclusive_visibility = is_visibility_exclusive;
}


void bdDataObject::SetDescription(const char *description)
{
    this->m_description.Clear();
    this->m_description.Append(description);
}


bdString& bdDataObject::GetDescription()
{
    return m_description;
}


int bdDataObject::IsLocked()
{
	return m_is_locked;
}

	
void bdDataObject::SetLocked(int is_locked)
{
	m_is_locked = is_locked;
}


bdDataObjectTag* bdDataObject::AddTag(bdString tag_namespace, bdString tag_class, bdString tag_object)
{
    bdDataObjectTag *tag = this->GetTag(tag_namespace, tag_class, tag_object);
    if(tag) return tag;
    
    tag = this->m_tags.AddNewToRightEnd();
    tag->m_namespace = tag_namespace;
    tag->m_class = tag_class;
    tag->m_object = tag_object;
    
    return tag;
}


bdDataObjectTag* bdDataObject::AddTag(const char *tag_namespace, const char *tag_class, const char *tag_object)
{
    bdString ns; ns.Assign(tag_namespace);
    bdString cl; cl.Assign(tag_class);
    bdString ob; ob.Assign(tag_object);
    return this->AddTag(ns, cl, ob);
}


bdDataObjectTag* bdDataObject::GetTag(bdString tag_namespace, bdString tag_class, bdString tag_object)
{
    bdListIterator<bdDataObjectTag> it;
    for(it.SetLeftEnd(this->m_tags); it.IsValid(); it.MoveRight())
    {
        bdDataObjectTag *tag = it.GetElementPointer();
        if((tag->m_namespace==tag_namespace) && (tag->m_class==tag_class) && (tag->m_object==tag_object)) { return tag; }
    }
    return NULL;
}


bdDataObjectTag* bdDataObject::GetTag(const char *tag_namespace, const char *tag_class, const char *tag_object)
{
    bdString ns; ns.Assign(tag_namespace);
    bdString cl; cl.Assign(tag_class);
    bdString ob; ob.Assign(tag_object);
    return this->GetTag(ns, cl, ob);
}


bdDataObjectTag* bdDataObject::GetTag()
{
    if(this->m_tags.IsEmpty()) return NULL;
    return this->m_tags.GetLeftEndPointerToElement();
}


int bdDataObject::AppendTagFrom(bdDataObject *obj)
{
    if(!obj) return 0;
    
    bdDataObjectTag *tag = this->GetTag();
    if(!tag) return 0;
    
    bdDataObjectTag *tag2 = obj->GetTag();
    if(!tag2) return 0;

    bdListIterator<bdString> iti, itv;
    for(iti.SetLeftEnd(tag2->GetListOfIdentifiers()), itv.SetLeftEnd(tag2->GetListOfValues()); iti.IsValid() && itv.IsValid(); iti.MoveRight(), itv.MoveRight())
    {
        tag->AddIdentifierAndValue(iti.GetElement(), itv.GetElement());
    }
    
    return 1;
}


int bdDataObject::RemoveTag(bdString tag_namespace, bdString tag_class, bdString tag_object)
{
    bdListNode<bdDataObjectTag> *tag_node_to_remove = NULL;
    bdListIterator<bdDataObjectTag> it;
    for(it.SetLeftEnd(this->m_tags); it.IsValid(); it.MoveRight())
    {
        bdDataObjectTag *tag = it.GetElementPointer();
        if((tag->m_namespace==tag_namespace) && (tag->m_class==tag_class) && (tag->m_object==tag_object))
        {
            tag_node_to_remove = it.GetNodeAddress();
            break;
        }
    }
    
    if(tag_node_to_remove)
    {
        this->m_tags.DeleteNode(tag_node_to_remove);
        return 1;
    }
    
    return 0;

}


int bdDataObject::RemoveTag(const char *tag_namespace, const char *tag_class, const char *tag_object)
{
    bdString ns; ns.Assign(tag_namespace);
    bdString cl; cl.Assign(tag_class);
    bdString ob; ob.Assign(tag_object);
    return this->RemoveTag(ns, cl, ob);
}


void bdDataObject::SaveTags(bdString &output_tags_string)
{
    output_tags_string.Clear();
    
    bdListIterator<bdDataObjectTag> it;
    for(it.SetLeftEnd(this->m_tags); it.IsValid(); it.MoveRight())
    {
        output_tags_string.Append("~");
        
        bdDataObjectTag *tag = it.GetElementPointer();
        output_tags_string.Append("<"); output_tags_string.Append(tag->m_namespace); output_tags_string.Append(">");
        output_tags_string.Append("["); output_tags_string.Append(tag->m_class); output_tags_string.Append("]");
        output_tags_string.Append("{"); output_tags_string.Append(tag->m_object); output_tags_string.Append("}");
        
        output_tags_string.Append("(");
        bdListIterator<bdString> its1; its1.SetLeftEnd(tag->m_identifiers);
        bdListIterator<bdString> its2; its2.SetLeftEnd(tag->m_values);
        while( its1.IsValid() && its2.IsValid() )
        {
            output_tags_string.Append("<"); output_tags_string.Append(its1.GetElement()); output_tags_string.Append(">");
            output_tags_string.Append("["); output_tags_string.Append(its2.GetElement()); output_tags_string.Append("]");
            its1.MoveRight();
            its2.MoveRight();
        }
        output_tags_string.Append(")");
        
        output_tags_string.Append("~");
    }
}


int bdDataObject::LoadTags(bdString &string_with_tags)
{
    this->m_tags.Reset();
    
    bdList<bdString> list_of_individual_tags;
    if(!string_with_tags.ExtractOptions('~', list_of_individual_tags)) return 0;
    if(list_of_individual_tags.IsEmpty()) return 0;
    
    bdListIterator<bdString> it_tags;
    for(it_tags.SetLeftEnd(list_of_individual_tags); it_tags.IsValid(); it_tags.MoveRight())
    {
        bdString ns, cl, ob, pa;
        
        if(!it_tags.GetElementPointer()->ExtractStringBetweenCharacters('<', '>', ns)) return 0;
        if(!it_tags.GetElementPointer()->ExtractStringBetweenCharacters('[', ']', cl)) return 0;
        if(!it_tags.GetElementPointer()->ExtractStringBetweenCharacters('{', '}', ob)) return 0;
        if(!it_tags.GetElementPointer()->ExtractStringBetweenCharacters('(', ')', pa)) return 0;
        
        bdDataObjectTag *tag = this->AddTag(ns, cl, ob);
        
        //cout<<"pa="<<pa.C_String()<<endl;
        
        bdList<bdString> list_of_identifiers;
        if(!pa.ExtractOptions('<', '>', list_of_identifiers)) return 0;
        if(list_of_identifiers.IsEmpty()) return 0;

        //cout<<"id_noel="<<list_of_identifiers.GetNumberOfElements()<<endl;
        
        bdList<bdString> list_of_values;
        if(!pa.ExtractOptions('[', ']', list_of_values)) return 0;
        if(list_of_values.GetNumberOfElements() != list_of_identifiers.GetNumberOfElements()) return 0;
        
        bdListIterator<bdString> it_ids;
        bdListIterator<bdString> it_values;
        it_ids.SetLeftEnd(list_of_identifiers);
        it_values.SetLeftEnd(list_of_values);
        while(it_ids.IsValid() && it_values.IsValid())
        {
            tag->AddIdentifierAndValue(it_ids.GetElement(), it_values.GetElement());
            
            it_ids.MoveRight();
            it_values.MoveRight();
        }
    }
    
    return 1;
}


int bdDataObject::LoadTags(const char *string_with_tags)
{
    bdString bds; bds.Assign(string_with_tags);
    return this->LoadTags(bds);
}


int bdDataObject::CopyFrom(bdDataObject *obj)
{
    if(!obj) return 0;
    
    this->m_history = obj->m_history;

    this->m_has_exclusive_visibility = obj->m_has_exclusive_visibility;

    this->m_description = obj->m_description;

    this->m_tags = obj->m_tags;
    
    // Do NOT copy object name, status (locked)
    //this->m_object_name = obj->m_object_name;
    //this->m_is_locked = obj->m_is_locked;

    this->m_color_RGBA[0] = obj->m_color_RGBA[0];
    this->m_color_RGBA[1] = obj->m_color_RGBA[1];
    this->m_color_RGBA[2] = obj->m_color_RGBA[2];
    this->m_color_RGBA[3] = obj->m_color_RGBA[3];

    return 1;
}



//-----------------------------------------------------------------------------------------------------------------------




bdFunctionObject::bdFunctionObject()
{
	this->SetupFunctionObject();
}


bdFunctionObject::~bdFunctionObject()
{
	this->Reset();
}


void bdFunctionObject::SetupFunctionObject()
{
	bdObject::SetupObject();
	m_is_abort_requested = 0;
}


void bdFunctionObject::SetAbortRequested(int is_abort_requested)
{
	m_is_abort_requested = is_abort_requested;
}


int bdFunctionObject::IsAbortRequested()
{
	return m_is_abort_requested;
}




//-----------------------------------------------------------------------------------------------------------------------





bdProcedureObject::bdProcedureObject()
{
    //this->SetupObject();
    //m_is_locked = 0;
    //this->m_color_RGBA[3] = 0;//set to no color.
    //this->m_has_exclusive_visibility = 0;// set to view multiple object simultaniously.
    //this->GetObjectName().Assign(" ");
    //this->AddTag("bd", "bdDataObject", " ");// basic tag that is always present when the object is created.
}


bdProcedureObject::~bdProcedureObject()
{
    
    //this->Reset();
    //m_is_locked = 0;
    //m_history.Clear();
}


//bdList<bdString>* bdProcedureObject::GetListOfRequestedDataObjectsByName()
//{
//    return &m_requested_data_objects_by_names;
//}


void bdProcedureObject::SetMain_3D_RenderWindowInteractor(vtkRenderWindowInteractor *main_3D_render_window_interactor)
{
    m_main_3D_render_window_interactor = main_3D_render_window_interactor;
}


int bdProcedureObject::AddRenderWindow(qbdVTKWidget *rw)
{
    if(!rw) return 0;
    this->m_render_windows.AddToRightEnd(rw);
    return 1;
}


bdList<qbdVTKWidget*>* bdProcedureObject::GetRenderWindows()
{
    return &m_render_windows;
}


int bdProcedureObject::SetReceivedDataObjects(bdList<bdDataObject*> &list_of_data_objects)
{
    if(list_of_data_objects.IsEmpty()) return 0;
    this->m_received_data_objects = list_of_data_objects;
    return 1;
}


bdList<bdDataObject*>* bdProcedureObject::GetReceivedDataObjects()
{
    return &(this->m_received_data_objects);
}


void bdProcedureObject::ResetListOfReceivedDataObjects()
{
    while(!this->m_received_data_objects.IsEmpty())
    {
        bdDataObject *dobj = this->m_received_data_objects.GetLeftEnd();
        dobj->SetLocked(0);
        this->m_received_data_objects.DeleteLeftEnd();
    }
}








//class bdProgressCounter
//{
////private:
//protected:
//
//	int m_value;//Current progress value.
//	int m_minimum;
//	int m_maximum;
//
////	bdObject *m_object;// Pointer to the object that updates this progress counter.
//
//public:
//
//	// Constructor / Destructor.
//	bdProgressCounter();
//	~bdProgressCounter();
//
////	// Set the object that will update the progress counter.
////	virtual void SetObject(bdObject *object);
//
////	// 
////	virtual bdObject* GetObject();
//
//
//	virtual void Set(int minimum, int maximum, int value);
//	virtual void Set(int minimum, int maximum){Set(minimum,maximum,minimum);};
//	virtual void SetMinimum(int min);
//	virtual void SetMaximum(int max);
//	virtual void SetValue(int value);
//
//	// Sets the value of the progress rescaled to the range max-min, where 'relative_value' is a number in range [0,100].
//	// The true value is calculated as true_value = min + (relative_value/100)*(max-min).
//	virtual void SetValueRelative(int relative_value);
//	
//	int GetMinimum();
//	int GetMaximum();
//	int GetValue();
//};
