#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_18_API_VERSION
#include <stdlib.h>
#include "Python.h"
#include "numpy/arrayobject.h"

int64_t bind_c_solve_2d_nonlinearconv_pyccel(int64_t n0_u, int64_t n1_u, double *u, int64_t n0_un, int64_t n1_un, double *un, int64_t n0_v, int64_t n1_v, double *v, int64_t n0_vn, int64_t n1_vn, double *vn, int64_t nt, double dt, double dx, double dy, int64_t c);

/*........................................*/



/*........................................*/

/*........................................*/
PyObject *solve_2d_nonlinearconv_pyccel_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    double *u;
    double *un;
    double *v;
    double *vn;
    int64_t nt;
    double dt;
    double dx;
    double dy;
    int64_t c;
    int64_t Out_0001;
    PyArrayObject *u_tmp;
    PyArrayObject *un_tmp;
    PyArrayObject *v_tmp;
    PyArrayObject *vn_tmp;
    PyObject *nt_tmp;
    PyObject *dt_tmp;
    PyObject *dx_tmp;
    PyObject *dy_tmp;
    PyObject *c_tmp;
    int64_t u_dim;
    int64_t u_dim_0001;
    int64_t un_dim;
    int64_t un_dim_0001;
    int64_t v_dim;
    int64_t v_dim_0001;
    int64_t vn_dim;
    int64_t vn_dim_0001;
    PyObject *result;
    static char *kwlist[] = {
        "u",
        "un",
        "v",
        "vn",
        "nt",
        "dt",
        "dx",
        "dy",
        "c",
        NULL
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!OOOOO", kwlist, &PyArray_Type, &u_tmp, &PyArray_Type, &un_tmp, &PyArray_Type, &v_tmp, &PyArray_Type, &vn_tmp, &nt_tmp, &dt_tmp, &dx_tmp, &dy_tmp, &c_tmp))
    {
        return NULL;
    }
    if (PyArray_NDIM(u_tmp) != 2)
    {
        PyErr_SetString(PyExc_TypeError, "\"u_tmp must have rank 2\"");
        return NULL;
    }
    else if (PyArray_TYPE(u_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"u must be double\"");
        return NULL;
    }
    else if (!PyArray_CHKFLAGS(u_tmp, NPY_ARRAY_C_CONTIGUOUS))
    {
        PyErr_SetString(PyExc_NotImplementedError, "\"Argument does not have the expected ordering (C)\"");
        return NULL;
    }
    else
    {
        u = PyArray_DATA(u_tmp);
    }
    if (PyArray_NDIM(un_tmp) != 2)
    {
        PyErr_SetString(PyExc_TypeError, "\"un_tmp must have rank 2\"");
        return NULL;
    }
    else if (PyArray_TYPE(un_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"un must be double\"");
        return NULL;
    }
    else if (!PyArray_CHKFLAGS(un_tmp, NPY_ARRAY_C_CONTIGUOUS))
    {
        PyErr_SetString(PyExc_NotImplementedError, "\"Argument does not have the expected ordering (C)\"");
        return NULL;
    }
    else
    {
        un = PyArray_DATA(un_tmp);
    }
    if (PyArray_NDIM(v_tmp) != 2)
    {
        PyErr_SetString(PyExc_TypeError, "\"v_tmp must have rank 2\"");
        return NULL;
    }
    else if (PyArray_TYPE(v_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"v must be double\"");
        return NULL;
    }
    else if (!PyArray_CHKFLAGS(v_tmp, NPY_ARRAY_C_CONTIGUOUS))
    {
        PyErr_SetString(PyExc_NotImplementedError, "\"Argument does not have the expected ordering (C)\"");
        return NULL;
    }
    else
    {
        v = PyArray_DATA(v_tmp);
    }
    if (PyArray_NDIM(vn_tmp) != 2)
    {
        PyErr_SetString(PyExc_TypeError, "\"vn_tmp must have rank 2\"");
        return NULL;
    }
    else if (PyArray_TYPE(vn_tmp) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_TypeError, "\"vn must be double\"");
        return NULL;
    }
    else if (!PyArray_CHKFLAGS(vn_tmp, NPY_ARRAY_C_CONTIGUOUS))
    {
        PyErr_SetString(PyExc_NotImplementedError, "\"Argument does not have the expected ordering (C)\"");
        return NULL;
    }
    else
    {
        vn = PyArray_DATA(vn_tmp);
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
    if (PyArray_IsScalar(dy_tmp, Float64))
    {
        PyArray_ScalarAsCtype(dy_tmp, &dy);
    }
    else if (PyFloat_Check(dy_tmp))
    {
        dy = PyFloat_AsDouble(dy_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"dy must be 64 bit real\"");
        return NULL;
    }
    if (PyArray_IsScalar(c_tmp, Int64))
    {
        PyArray_ScalarAsCtype(c_tmp, &c);
    }
    else if (PyLong_Check(c_tmp))
    {
        c = PyLong_AsLong(c_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"c must be 64 bit int\"");
        return NULL;
    }
    u_dim = PyArray_DIM(u_tmp, 0);
    u_dim_0001 = PyArray_DIM(u_tmp, 1);
    un_dim = PyArray_DIM(un_tmp, 0);
    un_dim_0001 = PyArray_DIM(un_tmp, 1);
    v_dim = PyArray_DIM(v_tmp, 0);
    v_dim_0001 = PyArray_DIM(v_tmp, 1);
    vn_dim = PyArray_DIM(vn_tmp, 0);
    vn_dim_0001 = PyArray_DIM(vn_tmp, 1);
    Out_0001 = bind_c_solve_2d_nonlinearconv_pyccel(u_dim, u_dim_0001, u, un_dim, un_dim_0001, un, v_dim, v_dim_0001, v, vn_dim, vn_dim_0001, vn, nt, dt, dx, dy, c);
    result = Py_BuildValue("l", Out_0001);
    return result;
}
/*........................................*/

static PyMethodDef mod_3hhsrbc7_methods[] = {
    {
        "solve_2d_nonlinearconv_pyccel",
        (PyCFunction)solve_2d_nonlinearconv_pyccel_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef mod_3hhsrbc7_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "mod_3hhsrbc7",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    -1,
    mod_3hhsrbc7_methods
};

/*........................................*/

PyMODINIT_FUNC PyInit_mod_3hhsrbc7(void)
{
    PyObject *m;

    import_array();

    m = PyModule_Create(&mod_3hhsrbc7_module);
    if (m == NULL) return NULL;

    return m;
}
