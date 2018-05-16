/**************************************************************************
 Object class as parent of all data analysis classes later on.
 
 Progress Counter class used in Function Object and Object classes to
 track the progress of method execution.
 
 Author: Danilo Babin
 File name: "bdObject.h"
***************************************************************************/



//To build as DLL, add:" /D "BD_OBJECT_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_OBJECT_EXPORTS
		#define BD_OBJECT_API __declspec(dllexport) 
	#else
		#define BD_OBJECT_API __declspec(dllimport) 
	#endif
#else // GCC
    #ifdef BD_OBJECT_EXPORTS
        #define BD_OBJECT_API __attribute__((visibility("default")))
	#else
		#define BD_OBJECT_API
	#endif
#endif




#ifndef BD_OBJECT_DEF
	#define BD_OBJECT_DEF


#include "bdString.h"
#include <vector>



class bdObject;


class BD_OBJECT_API bdProgressCounter
{
protected:
    
    /// Current progress value.
    int m_value;
    
    /// Minimum progress value (start value).
    int m_minimum;
    
    /// Maximum progress value (end value).
    int m_maximum;
    
    /// Pointer to the object that updates this progress counter.
    bdObject *m_object;
    
    friend class bdObject;
    
    /// Method that if a FRIEND to bdObject and sets its pointer to bdProgressCounter to NULL.
    /// Called in destructors to break the connection between objects.
    void DisconnectObject();
    
public:
    
    /// Constructor / Destructor.
    bdProgressCounter();
    ~bdProgressCounter();
    
    /// Get the object that updates the progress counter.
    bdObject* GetObject();
    
    /// Set the range (minimum and maximum) and value of the progress counter.
    virtual void Set(int minimum, int maximum, int value);
    virtual void Set(int minimum, int maximum){Set(minimum,maximum,minimum);};
    virtual void SetMinimum(int min);
    virtual void SetMaximum(int max);
    virtual void SetValue(int value);
    
    /// Sets the value of the progress rescaled to the range max-min, where 'relative_value' is a number in range [0,100].
    /// The true value is calculated as true_value = min + (relative_value/100)*(max-min).
    virtual void SetValueRelative(int relative_value);
    
    /// Get range min and max.
    int GetMinimum();
    int GetMaximum();
    
    /// Get value of the counter.
    int GetValue();
};






//--------------------------------------------------------------------------------------------------




class BD_OBJECT_API bdObject
{
private:

	/// Progress counter object that tracks the progress of method execution.
	bdProgressCounter *m_progress_counter;

	friend void bdProgressCounter::DisconnectObject();

protected:

	/// Method that does the work of the Object constructor. It is made so that it could be 
	/// called from constructors of inheriting objects.
	virtual void SetupObject();
    
    /// String containing the type name.
    bdString m_type_name;

public:

	/// Constructor.
	bdObject();

	/// Destructor.
	~bdObject();
    
    /// Get string containing the type name.
    bdString& GetTypeName();

	/// Reset the object.
	virtual void Reset();

	/// Connect a progress counter object to track execution progress.
	int ConnectProgressCounter(bdProgressCounter *progress_counter);

	/// Method that is a FRIEND to bdProgressCounter and sets its pointer to Object to NULL.
	/// Called in destructors to break the connection between objects.
	int DisconnectProgressCounter();

	/// Set properties of a progress counter.
	virtual void SetProgressCounterProperties(int min, int max, int value);

	/// Calls SetValueRelative() member function of bdProgressCounter. 
	virtual void SetProgressCounterRelativeValue(int relative_value);

};



//----------------------------------------------------------------------------------------------------


// Pre-declare some classes for later use
class vtkRenderWindowInteractor;
class vtkActor;
class qbdDataObjectContextMenu;
class qbdDataObjectPropertiesWidget;
class qbdDataContainer;
class bdDataObject;
class QwtPlot;
class QwtPlotZoomer;



typedef void (bdDataObject::*bdDataObjectMethod)(qbdDataContainer *data_container);


/// Helper template class of batch call mapper node for the list of mapper commands.

template<class T>
class bdObjectBatchCallMapperCommandNode
{
public:
    bdString m_command;
    bdString m_parameters;
    bdString m_help_text;
    T m_method_pointer; // bdDataObjectMethod m_method;
};


/// Helper template class of batch call mapper.

template<class T>
class bdObjectBatchCallMapper
{
    bdList< bdObjectBatchCallMapperCommandNode<T> > m_batch_command_map;
    
public:
    
    /// Search the m_batch_command_map and for a given command return the pointer to the method assigned to it.
    T GetMethod(const char *command);
    
    /// Make a new entry in the m_batch_command_map.
    int RegisterMethod(const char *command, const char *parameters, const char *help_text, T method_pointer);
};

template<class T>
T bdObjectBatchCallMapper<T>::GetMethod(const char *command)
{
    bdString str_command; str_command.Append(command);
    bdListIterator< bdObjectBatchCallMapperCommandNode<T> > it;
    for(it.SetLeftEnd(m_batch_command_map); it.IsValid(); it.MoveRight())
    {
        if(it.GetElementPointer()->m_command == str_command)
        {
            return it.GetElementPointer()->m_method_pointer;
        }
    }
    return NULL;
};

template<class T>
int bdObjectBatchCallMapper<T>::RegisterMethod(const char *command, const char *parameters, const char *help_text, T method_pointer)
{
    bdObjectBatchCallMapperCommandNode<T> *new_element = m_batch_command_map.AddNewToRightEnd();
    new_element->m_command.Append(command);
    new_element->m_parameters.Append(parameters);
    new_element->m_help_text.Append(help_text);
    new_element->m_method_pointer = method_pointer;
    return 1;
};



/// Helper class of data object tags.

class BD_OBJECT_API bdDataObjectTag
{
public:
    bdString m_namespace;
    bdString m_class;
    bdString m_object;
    bdList<bdString> m_identifiers;
    bdList<bdString> m_values;
    
    
    bdList<bdString>* GetListOfValues(){return &m_values;};
    bdList<bdString>* GetListOfIdentifiers(){return &m_identifiers;};
    
    bdString GetValueForGivenIdentifier(bdString identifier);
    bdString GetValueForGivenIdentifier(const char *identifier);

    int GetValueAsIntegerForGivenIdentifier(bdString identifier, int &output_integer);
    int GetValueAsIntegerForGivenIdentifier(const char *identifier, int &output_integer);

    int GetValueAsDoubleForGivenIdentifier(bdString identifier, double &output_double);
    int GetValueAsDoubleForGivenIdentifier(const char *identifier, double &output_double);

    void AddIdentifierAndValue(bdString identifier, bdString value);
    void AddIdentifierAndValue(const char *identifier, const char *value);

    /// Add reference to data object given by 'data_object_name' (the string will have '*' added in front) and value for sub-object.
    void AddDataObjectReference(bdString data_object_name, bdString value);
    void AddDataObjectReference(const char *data_object_name, const char *value);
    
    /// Get value for given data object name.
    int GetValueForGivenDataObjectName(bdString data_object_name, bdString &output_value);
    int GetValueForGivenDataObjectName(const char *data_object_name, bdString &output_value);
    
    /// Check if the refence to a certain object and a value exists.
    int HasDataObjectReference(bdString data_object_name, bdString value);
    int HasDataObjectReference(const char *data_object_name, const char *value);

    /// Check if the refence to a certain object exists (no tag value is used in search, only identifier).
    int HasDataObjectReference(bdString data_object_name);
    int HasDataObjectReference(const char *data_object_name);
    
    
    
    bdDataObjectTag& operator =(bdDataObjectTag &tag);    
};



/// Base class for data objects.

class BD_OBJECT_API bdDataObject : public bdObject
{
protected:

	/// String containing signatures(names) and parameters of all functions called on this object.
	bdString m_history;

	/// Indicator showing if the object is locked (used) by an algorithm, so it should not be used for another function at the
	/// same time. Notice: This indicator does NOT stop a method/function being executed on the object, it should only be set by
	/// a function to indicate that it is being used and should be reset after the function is completed.
	/// This indicator will be used on Qt level to allow or refuse the appearance of context menu.
	int m_is_locked;

    /// Indicator showing if the object canbe viewed together with other objects of this type. If this value is non-zero, only a single object can be viewed at a time.
    int m_has_exclusive_visibility;
    
//    /// String containing the type name.
//    bdString m_type_name;
    
    /// String containing short description.
    bdString m_description;
    
    /// Tag is a string added by algorithms to indicate that they should be used as input to another method (to communicate
    /// data between functions on the GUI level).
    //bdString m_tag;
    bdList<bdDataObjectTag> m_tags;

    /// List of comments inserted by various algorithms. They are of type:
    /// <namespace>[class]{object}(<identifier1>[value1]<identifier2>[value2]...)
    /// example: <bd>[PWV]{aortic_position}(<coordinates>[13,3,78.2])
    //bdList<bdString> m_algorithm_comments;
    
//    /// In case an object contains multiple sub-objects, a tag can be assinged to each subobject by entering a subtag to this list of subtags.
//    bdList<bdString> m_subtags;
    
    /// String containing object name.
    bdString m_object_name;
    
    /// Color of the data object (to be used optionally if needed).
    unsigned char m_color_RGBA[4];
    
    
    /// NOTE: Child class (the one that inherits for this object) will have a pointer to its context menu (that will inherit from qbdDataObjectContextMenu in "qbdObject.h") and
    /// a pointer to its properties widget (that will inherit from qbdDataObjectPropertiesWidget in "qbdObject.h").
    
    /// NOTE: Child class (the one that inherits for this object) will have a map between batch commands and actual methods of the class. In this way we will be able to call
    /// the methods of each class through batch command line. Declaration of mapper list will be in the following fashion:
    /// bdObjectBatchCallMapper<TypeName_Of_Pointer_To_Method>  m_batch_call_mapper;

public:

	/// Constructor.
	bdDataObject();

	/// Destructor.
	~bdDataObject();

	/// Get history string (to track called functions on this object).
	bdString& GetHistory();
    
//    /// Get string containing the type name.
//    bdString& GetTypeName();
    
    /// Get string containing the object name.
    bdString& GetObjectName();
    
    /// Check if multiple objects can be viewed simultaneously.
    int HasExclusiveVisibility();

    /// Set exclusive visibility.
    void SetExclusiveVisibility(int is_visibility_exclusive);

    /// Set the description string.
    void SetDescription(const char *description);
    
    /// Get the object color.
    unsigned char GetColor_Red(){return m_color_RGBA[0];};
    unsigned char GetColor_Green(){return m_color_RGBA[1];};
    unsigned char GetColor_Blue(){return m_color_RGBA[2];};
    unsigned char GetColor_A(){return m_color_RGBA[3];};
    
    /// Get the description string.
    bdString& GetDescription();

	/// Check if the object is being used by a function (so it should not be used by another function/method at the same time).
	int IsLocked();

	/// Lock/unlock the object. Use this to indicate that an external function is using the object. Do not forget to release the 
	/// object after the efunction execution has ended.
	void SetLocked(int is_locked);
    
    
    /// Add new tag with given tag namespace, class and object name and return pointer to it. If tag already exists just returns pointer to it.
    bdDataObjectTag* AddTag(bdString tag_namespace, bdString tag_class, bdString tag_object);
    bdDataObjectTag* AddTag(const char *tag_namespace, const char *tag_class, const char *tag_object);

    /// Get tag with given tag namespace, class and object name. If none found, return NULL.
    bdDataObjectTag* GetTag(bdString tag_namespace, bdString tag_class, bdString tag_object);
    bdDataObjectTag* GetTag(const char *tag_namespace, const char *tag_class, const char *tag_object);
    
    /// Get the first tag in the list of tags.
    bdDataObjectTag* GetTag();
    
    /// Copy tag values of the default (first) tag from input object.
    int AppendTagFrom(bdDataObject *obj);

    
    /// Remove tag. If it does not exist return fail 0, success 1.
    int RemoveTag(bdString tag_namespace, bdString tag_class, bdString tag_object);
    int RemoveTag(const char *tag_namespace, const char *tag_class, const char *tag_object);
    
    /// Save tags of the object. Tag should be of following type:
    /// ~<namespace>[class]{object}(<identifier1>[value1]<identifier2>[value2]...)~
    /// example: ~<bd>[PWV]{velocity_curve}(<index>[0]<version>[1])~
    void SaveTags(bdString &output_tags_string);
    int LoadTags(const char *string_with_tags);
    int LoadTags(bdString &string_with_tags);

    
    
    /// Set the ploter and zoomer for the 2D signal plot.
    virtual void SetPlotWidgetAndZoomer(QwtPlot *plot_widget, QwtPlotZoomer *zoomer){cout<<(m_type_name.C_String())<<" called bdDataObject::SetPlotWidgetAndZoomer()"<<endl;};
    
    /// Set the interactor for VTK visualization.
    virtual void SetInteractor(vtkRenderWindowInteractor *interactor){cout<<(m_type_name.C_String())<<" called bdDataObject::SetInteractor()"<<endl;};

    /// Get the interactor for VTK visualization.
    virtual vtkRenderWindowInteractor* GetInteractor(){cout<<(m_type_name.C_String())<<" called bdDataObject::GetInteractor()"<<endl; return NULL;};
    
    /// Set the visibility of the object for VTK visualization.
    virtual void SetVisibility(int is_visible){cout<<(m_type_name.C_String())<<" called bdDataObject::SetVisibility()"<<endl;};
    
    /// Check if the actor belongs to the object (used for calling the context menu).
    virtual int HasActor(vtkActor *actor){cout<<(m_type_name.C_String())<<" called bdDataObject::HasActor()"<<endl; return 0;};
    
    /// Get the context menu of the object.
    virtual qbdDataObjectContextMenu* GetContextMenu(){cout<<(m_type_name.C_String())<<" called bdDataObject::GetContextMenu()"<<endl; return NULL;};
    
    /// Get the properties widget of the object.
    virtual qbdDataObjectPropertiesWidget* GetPropertiesWidget(){cout<<(m_type_name.C_String())<<" called bdDataObject::GetPropertiesWidget()"<<endl; return NULL;};
    
    /// Update the visual structures of object after modification.
    virtual void Update(){cout<<(m_type_name.C_String())<<" called bdDataObject::Update()"<<endl;};
    
    /// Method that will do the saving to a default file type.
    virtual void SaveToDefaultFileType(const char *file_path){cout<<(m_type_name.C_String())<<" called bdDataObject::SaveToDefaultFileType()"<<endl;};
    
    /// Make an exact copy from input object.
    virtual int CopyFrom(bdDataObject *obj);//{cout<<(m_type_name.C_String())<<" called bdDataObject::CopyFrom()"<<endl; return 0;};
    
    /// This will show context menu of the object.
    virtual void ShowContextMenu(){cout<<(m_type_name.C_String())<<" called bdDataObject::ShowContextMenu()"<<endl;};

    /// Method that will show basic context menu of the object.
    virtual void ShowBasicContextMenu(){cout<<(m_type_name.C_String())<<" called bdDataObject::ShowBasicContextMenu()"<<endl;};
    
    /// Method that will show the object properties in the properties widget.
    virtual void ShowProperties(){cout<<(m_type_name.C_String())<<" called bdDataObject::ShowProperties()"<<endl;};
    
    /// Method that will execute batch command given by input string and data object.
    virtual int ProcessBatchCommand(const char *command, qbdDataContainer *data_container){cout<<(m_type_name.C_String())<<" called bdDataObject::ProcessBatchCommand()"<<endl; return 0;};
    
    /// This will calculate the bounds of the object in 3D space. If it's a 2D object (e.g. a signal) the return value should be 0 (3D bounds not applicable).
    virtual int Bounds_3D_WorldCoordinates(double &x_min, double &x_max, double &y_min, double &y_max, double &z_min, double &z_max)
    {x_min = y_min = z_min = 0; x_max = y_max = z_max = 100; cout<<(m_type_name.C_String())<<" called bdDataObject::Bounds_3D_WorldCoordinates()"<<endl; return 0;};
    
    /// Check if the object is visible.
    virtual int IsVisible(){cout<<(m_type_name.C_String())<<" called bdDataObject::IsVisible()"<<endl; return 0;};

    
};



//LATEST

template<class T>
int Convert(bdDataObject *data_object, T* &out_child_object)
{
    T child;
    if(data_object->GetTypeName() == child.GetTypeName())
    {
        out_child_object = dynamic_cast<T*>(data_object);
        //T *child_object = dynamic_cast<T*>(data_object);
        return 1;
    }
    return 0;
};


template<class T>
int Convert(std::vector<bdDataObject*> data_objects, std::vector<T*> &out_child_objects)
{
    T child;
    int return_value = 0;
    for(int i=0; i<data_objects.size(); i++)
    {
        //cout<<(data_objects[i]->GetTypeName())<<"_"<<child.GetTypeName()<<endl;
        if(data_objects[i]->GetTypeName() == child.GetTypeName())
        {
            //cout<<" IN! "<<endl;
            T *child_object = dynamic_cast<T*>(data_objects[i]);
            out_child_objects.push_back(child_object);
            return_value = 1;
        }
    }
    return return_value;
};


template<class T1, class T2>
void Convert(std::vector<T1*> &vector_to_convert, bdList<T2*> &out_converted_list)
{
    out_converted_list.Reset();
    for(int i=0; i<vector_to_convert.size(); i++)
    {
        T2 *p = (T2*) vector_to_convert[i];
        out_converted_list.AddToRightEnd(p);
    }
};



//template<class T>
//int ConvertToParent(bdDataObject *data_object, T *out_child_object)
//{
//    T child;
//    if(data_object->GetTypeName() == child->GetTypeName())
//    {
//        T *child_object = dynamic_cast<T*>(data_object);
//        return 1;
//    }
//    return 0;
//};


//template<class T>
//int Convert(std::vector<bdDataObject*> data_objects, std::vector<T*> out_child_objects)
//{
////    if(data_objects.size()<=0) return 0;
//    T child;
//    for(int i=0; i<data_objects.size(); i++)
//    {
//        if(data_objects[i]->GetTypeName() == child->GetTypeName())
//        {
//            T *child_object = dynamic_cast<T*>(data_objects[i]);
//            return 1;
//        }
//    }
//    return 0;
//};


// Convert(data_object, crv);



//----------------------------------------------------------------------------------------------------




class BD_OBJECT_API bdFunctionObject : public bdObject
{
private:

	/// ABORT flag of the function object. If set indicates that the method should abort ASAP.
	int m_is_abort_requested;

protected:

	/// Initializes the object.
	virtual void SetupFunctionObject();


public:

	/// Constructor.
	bdFunctionObject();

	/// Destructor.
	~bdFunctionObject();

	/// Set the ABORT flag of the function object.
	void SetAbortRequested(int is_abort_requested);

	/// Chack if the ABORT flag is set.
	int IsAbortRequested();

};



//----------------------------------------------------------------------------------------------------



// Pre-declare needed classes.
class qbdProcedureObjectPropertiesWidget;
class qbdVTKWidget;
//class vtkRenderWindow;
//class vtkCamera;



class BD_OBJECT_API bdProcedureObject : public bdObject
{
private:
        
    /// String containing short description, usually object name.
    bdString m_description;
    
    /// Child class (the one that inherits for this object) will have a pointer to its context menu (that will inherit from qbdDataObjectContextMenu in "qbdObject.h") and
    /// a pointer to its properties widget (that will inherit from qbdDataObjectPropertiesWidget in "qbdObject.h").
    
    
protected:
    
    /// Data objects that are created during procedure execution which can be added as results to the application.
    bdList<bdDataObject> m_objects_for_export;
    
//    /// The list of data requested by the procedure. Normally, the procedure properties widget will contain a input table to pass data names. When the user enters a name into the table,
//    /// the list will be updated (this is done in the child of procedure properties widget).
//    bdList<bdString> m_requested_data_objects_by_names;
    
    /// Data objects that were received from the main app for execution of the procedure. The data objects are locked by
    /// the main app, so this list is also used to unlock required data after procedure has ended.
    bdList<bdDataObject*> m_received_data_objects;
    
    /// Requested render windows from the main app.
    bdList<qbdVTKWidget*> m_render_windows;
    
    /// The main 3D render window interactor.
    vtkRenderWindowInteractor *m_main_3D_render_window_interactor;
    
// SHOULD ALSO ADD POINTERS TO MAIN SIGNAL PLOT WIDGET AND LIST OF REQUESTED NEW SIGNAL PLOT WIDGETS (MAYBE ALSO FOR TABLES).
    
    
    
public:
    
    /// Constructor.
    bdProcedureObject();
    
    /// Destructor.
    ~bdProcedureObject();
    
    
    
    /// Get the number of requested render windows by the procedure. Normally, the procedure might require additional render windows that need to be created by the main app. Hence, when the procedure is created,
    /// the main app gets the number of render windows that it needs, creates these render windows and passes them to the procedure. When the procedure ends, the main app removes (deallocates) these render windows.
    virtual int GetNumberOfRequestedRenderWindows() { cout<<(m_type_name.C_String())<<" called bdProcedureObject::GetNumberOfRequestedRenderWindows()"<<endl; return 0;};

    /// Get the properties widget of the object.
    virtual qbdProcedureObjectPropertiesWidget* GetPropertiesWidget(){cout<<(m_type_name.C_String())<<" called bdProcedureObject::GetPropertiesWidget()"<<endl; return NULL;};
    
    /// Update the visual structures of procedure after modification (e.g. after new data has been added or selected...). Return fail 0, success 1. 
    virtual int Update(){cout<<(m_type_name.C_String())<<" called bdProcedureObject::Update()"<<endl; return 0;};
    
    /// Perform required data release at closing of the procedure.
    virtual void Close(){cout<<(m_type_name.C_String())<<" called bdProcedureObject::Close()"<<endl;};
    
    /// Get the list of data requested by the procedure. Normally, the procedure properties widget will contain a input table to pass data names. When the request to get the specified data is made (by the user),
    /// this method will return the list of data requested and the main app should respond to this by passing data pointers to requested data objects.
    virtual void GetListOfRequestedDataObjectsByName(bdList<bdString> *output_list_of_names){cout<<(m_type_name.C_String())<<" called bdProcedureObject::GetListOfRequestedDataObjectsByName()"<<endl;};
    
    
    
    /// Set the pointer to the main 3D render window.
    void SetMain_3D_RenderWindowInteractor(vtkRenderWindowInteractor *main_3D_render_window_interactor);
    
    /// Add an existing render window to this procedure (procedure makes a request for a render window, the main app created it and passes it to procedure using this method).
    int AddRenderWindow(qbdVTKWidget *rw);
    
    /// Get the list of added render windows (created by the main app and used by this procedure).
    bdList<qbdVTKWidget*>* GetRenderWindows();
    
    /// Set the data objects that were received by the main app for execution of the procedure. The data objects are locked by
    /// the main app, so this list is also used to unlock required data after procedure has ended.
    int SetReceivedDataObjects(bdList<bdDataObject*> &list_of_data_objects);
    
    /// Get the pointer to the list of received data objects.
    bdList<bdDataObject*>* GetReceivedDataObjects();
    
    /// Removes all data object pointers froom the list of received objects and unocks them.
    void ResetListOfReceivedDataObjects();
    
    
    /// Get the interactor for VTK visualization.
    //    virtual vtkRenderWindowInteractor* GetInteractor(){cout<<(m_type_name.C_String())<<" called bdDataObject::GetInteractor()"<<endl; return NULL;};
    
//    /// Set the visibility of the object for VTK visualization.
//    virtual void SetVisibility(int is_visible){cout<<(m_type_name.C_String())<<" called bdDataObject::SetVisibility()"<<endl;};
    
    
//    /// Method that will show the procedure properties in the properties widget.
//    virtual void ShowProperties(){cout<<(m_type_name.C_String())<<" called bdDataObject::ShowProperties()"<<endl;};

    
};


#endif