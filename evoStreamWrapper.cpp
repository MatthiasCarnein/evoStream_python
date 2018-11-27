#include <Python.h>
#include <Vector>

#include <cstdio>

#include "evoStream.hpp"


typedef struct {
    PyObject_HEAD
    EvoStream * evoStream;
} PyEvoStream;




static PyModuleDef evoStreammodule = {
    PyModuleDef_HEAD_INIT,
    "evoStream",
    "Example module that wrapped a C++ object",
    -1,
    NULL, NULL, NULL, NULL, NULL
};

// initialize PyEvoStream Object
static int PyEvoStream_init(PyEvoStream *self, PyObject *args, PyObject *kwds){
    double r;
    double lambda;
    int tgap;
    int k;
    double crossoverRate;
    double mutationRate;
    int populationSize;
    int initializeAfter;
    int reclusterGenerations;

    if (!PyArg_ParseTuple(args, "ddiiddiii", &r, &lambda, &tgap, &k, &crossoverRate, &mutationRate, &populationSize, &initializeAfter, &reclusterGenerations))
        return -1;

    self->evoStream = new EvoStream(r, lambda, tgap, k, crossoverRate, mutationRate, populationSize, initializeAfter, reclusterGenerations);

    return 0;
}

// destruct the object
static void PyEvoStream_dealloc(PyEvoStream * self){
    delete self->evoStream;
    Py_TYPE(self)->tp_free(self);
}


// cluster
static PyObject *PyEvoStream_cluster(PyEvoStream *self, PyObject *args)
{
    PyObject *float_list;
    int size;

    if (!PyArg_ParseTuple(args, "O", &float_list))
        return NULL;

    size = PyObject_Length(float_list);
    if (size < 0)
        return NULL;

    std::vector<double> data(size);

    for (int i = 0; i < size; i++) {
        PyObject *item;
        item = PyList_GetItem(float_list, i);
        data[i] = PyFloat_AsDouble(item);
    }
    self->evoStream->cluster(data);
    return Py_BuildValue("i", 0);
}



// get_microclusters
static PyObject *PyEvoStream_get_microclusters(PyEvoStream *self, PyObject *args){
    std::vector<std::vector<double> > clusters = self->evoStream->get_microclusters();

    PyObject *rows = PyList_New(0);
    for(int i=0; i < clusters.size(); i++){
        PyObject *cols = PyList_New(0);
        for(int j=0; j < clusters[0].size(); j++){
            PyList_Append(cols, PyFloat_FromDouble(clusters[i][j]));
        }
        PyList_Append(rows, cols);
    }
    return rows;
}


// get_microweights
static PyObject *PyEvoStream_get_microweights(PyEvoStream *self, PyObject *args){
    std::vector<double> weights = self->evoStream->get_microweights();

    PyObject *newlist = PyList_New(0);
    for(int i=0; i < weights.size(); i++){
        PyList_Append(newlist, PyFloat_FromDouble(weights[i]));
    }
    return newlist;
}


// get_macroclusters
static PyObject *PyEvoStream_get_macroclusters(PyEvoStream *self, PyObject *args){
    std::vector<std::vector<double> > clusters = self->evoStream->get_macroclusters();

    PyObject *rows = PyList_New(0);
    for(int i=0; i < clusters.size(); i++){
        PyObject *cols = PyList_New(0);
        for(int j=0; j < clusters[0].size(); j++){
            PyList_Append(cols, PyFloat_FromDouble(clusters[i][j]));
        }
        PyList_Append(rows, cols);
    }
    return rows;
}


// get_macroweights
static PyObject *PyEvoStream_get_macroweights(PyEvoStream *self, PyObject *args){
    std::vector<double> weights = self->evoStream->get_macroweights();

    PyObject *newlist = PyList_New(0);
    for(int i=0; i < weights.size(); i++){
        PyList_Append(newlist, PyFloat_FromDouble(weights[i]));
    }
    return newlist;
}

// microToMacro
static PyObject *PyEvoStream_microToMacro(PyEvoStream *self, PyObject *args){

    std::vector<int> assignment = self->evoStream->microToMacro();

    PyObject *newlist = PyList_New(0);
    for(int i=0; i < assignment.size(); i++){
        PyList_Append(newlist, PyLong_FromLong(assignment[i]));
    }
    return newlist;
}


// recluster
static PyObject *PyEvoStream_recluster(PyEvoStream *self, PyObject *args){

    int generations;
    PyArg_ParseTuple(args, "i", &generations);

    self->evoStream->recluster(generations);

    return Py_BuildValue("i", 0);
}







// methods
static PyMethodDef PyEvoStream_methods[] = {
    { "cluster", (PyCFunction)PyEvoStream_cluster,    METH_VARARGS,       "insert observation" },
    { "get_microclusters", (PyCFunction)PyEvoStream_get_microclusters,    METH_VARARGS,       "get centres of micro clusters" },
    { "get_microweights", (PyCFunction)PyEvoStream_get_microweights,    METH_VARARGS,       "get weight of micro clusters" },
    { "get_macroweights", (PyCFunction)PyEvoStream_get_macroweights,    METH_VARARGS,       "get weights of macro clusters" },
    { "get_macroclusters", (PyCFunction)PyEvoStream_get_macroclusters,    METH_VARARGS,       "get centres of macro clusters" },
    { "microToMacro", (PyCFunction)PyEvoStream_microToMacro,    METH_VARARGS,       "assignment of micro clusters to macro clusters" },
    { "recluster", (PyCFunction)PyEvoStream_recluster,    METH_VARARGS,       "start reclustering" },
    {NULL}  /* Sentinel */
};

static PyTypeObject PyEvoStreamType = { PyVarObject_HEAD_INIT(NULL, 0)
                                    "evoStream.EvoStream"   /* tp_name */
                                };


// create the module
PyMODINIT_FUNC PyInit_evoStream(void){
    PyObject* m;

    PyEvoStreamType.tp_new = PyType_GenericNew;
    PyEvoStreamType.tp_basicsize=sizeof(PyEvoStream);
    PyEvoStreamType.tp_dealloc=(destructor) PyEvoStream_dealloc;
    PyEvoStreamType.tp_flags=Py_TPFLAGS_DEFAULT;
    PyEvoStreamType.tp_doc="EvoStream objects";
    PyEvoStreamType.tp_methods=PyEvoStream_methods;
    //~ PyEvoStreamType.tp_members=Noddy_members;
    PyEvoStreamType.tp_init=(initproc)PyEvoStream_init;

    if (PyType_Ready(&PyEvoStreamType) < 0)
        return NULL;

    m = PyModule_Create(&evoStreammodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&PyEvoStreamType);
    PyModule_AddObject(m, "EvoStream", (PyObject *)&PyEvoStreamType); // Add EvoStream object to the module
    return m;
}
