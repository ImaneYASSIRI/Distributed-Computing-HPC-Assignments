#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_18_API_VERSION
#include <stdlib.h>
#include "Python.h"
#include "numpy/arrayobject.h"

int64_t bind_c_solve_1d_burger_pyccel(int64_t n0_u, double *u, int64_t n0_un, double *un, int64_t nt, int64_t nx, double dt, double dx, double nu);

/*........................................*/



/*........................................*/

/*........................................*/
PyObject *solve_1d_burger_pyccel_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double *u;
    double *un;
    int64_t nt;
    int64_t nx;
    double dt;
    double dx;
    double nu;
    int64_t Out_0001;
    PyArrayObject *u_tmp;
    PyArrayObject *un_tmp;
    PyObject *nt_tmp;
    PyObject *nx_tmp;
    PyObject *dt_tmp;
    PyObject *dx_tmp;
    PyObject *nu_tmp;
    int64_t u_dim;
    int64_t un_dim;
    PyObject *result;
    static char *kwlist[] = {
        "u",
        "un",
        "nt",
        "nx",
        "dt",
        "dx",
        "nu",
        NULL
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!OOOOO", kwlist, &PyArray_Type, &u_tmp, &PyArray_Type, &un_tmp, &nt_tmp, &nx_tmp, &dt_tmp, &dx_tmp, &nu_tmp))
    {
        return NULL;
    }
    if (PyArray_NDIM(u_tmp) != 1)
    {
        PyErr_SetString(PyExc_TypeError, "\"u_tmp must have rank 1\"");
        return NULL;
    }
    else if (PyArray_TYPE(u_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"u must be double\"");
        return NULL;
    }
    else
    {
        u = PyArray_DATA(u_tmp);
    }
    if (PyArray_NDIM(un_tmp) != 1)
    {
        PyErr_SetString(PyExc_TypeError, "\"un_tmp must have rank 1\"");
        return NULL;
    }
    else if (PyArray_TYPE(un_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"un must be double\"");
        return NULL;
    }
    else
    {
        un = PyArray_DATA(un_tmp);
    }
    if (PyArray_IsScalar(nt_tmp, Int64))
    {
        PyArray_ScalarAsCtype(nt_tmp, &nt);
    }
    else if (PyLong_Check(nt_tmp))
    {
        nt = PyLong_AsLong(nt_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"nt must be 64 bit int\"");
        return NULL;
    }
    if (PyArray_IsScalar(nx_tmp, Int64))
    {
        PyArray_ScalarAsCtype(nx_tmp, &nx);
    }
    else if (PyLong_Check(nx_tmp))
    {
        nx = PyLong_AsLong(nx_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"nx must be 64 bit int\"");
        return NULL;
    }
    if (PyArray_IsScalar(dt_tmp, Float64))
    {
        PyArray_ScalarAsCtype(dt_tmp, &dt);
    }
    else if (PyFloat_Check(dt_tmp))
    {
        dt = PyFloat_AsDouble(dt_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"dt must be 64 bit real\"");
        return NULL;
    }
    if (PyArray_IsScalar(dx_tmp, Float64))
    {
        PyArray_ScalarAsCtype(dx_tmp, &dx);
    }
    else if (PyFloat_Check(dx_tmp))
    {
        dx = PyFloat_AsDouble(dx_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"dx must be 64 bit real\"");
        return NULL;
    }
    if (PyArray_IsScalar(nu_tmp, Float64))
    {
        PyArray_ScalarAsCtype(nu_tmp, &nu);
    }
    else if (PyFloat_Check(nu_tmp))
    {
        nu = PyFloat_AsDouble(nu_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"nu must be 64 bit real\"");
        return NULL;
    }
    u_dim = PyArray_DIM(u_tmp, 0);
    un_dim = PyArray_DIM(un_tmp, 0);
    Out_0001 = bind_c_solve_1d_burger_pyccel(u_dim, u, un_dim, un, nt, nx, dt, dx, nu);
    result = Py_BuildValue("l", Out_0001);
    return result;
}
/*........................................*/

static PyMethodDef mod_lxhtwhvo_methods[] = {
    {
        "solve_1d_burger_pyccel",
        (PyCFunction)solve_1d_burger_pyccel_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef mod_lxhtwhvo_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "mod_lxhtwhvo",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    -1,
    mod_lxhtwhvo_methods
};

/*........................................*/

PyMODINIT_FUNC PyInit_mod_lxhtwhvo(void)
{
    PyObject *m;

    import_array();

    m = PyModule_Create(&mod_lxhtwhvo_module);
    if (m == NULL) return NULL;

    return m;
}
