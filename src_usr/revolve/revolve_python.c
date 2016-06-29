#include <python2.7/Python.h>

extern int wrap_revolve(int* check,int* capo,int* fine,int *snaps_in,int* info); 

static PyObject* revolve_whatodo(PyObject *self, PyObject *args)
{
  int check, capo, fine, snaps_in, info, whatodo;

  if (!PyArg_ParseTuple(args, "iiiii", &check, &capo, &fine, &snaps_in, &info))
    return NULL;
  whatodo = wrap_revolve(&check, &capo, &fine, &snaps_in, &info);
  return Py_BuildValue("iiiii", whatodo,check,capo,fine,info);
}

static PyMethodDef RevolveMethods[] = {
            {"whatodo",  revolve_whatodo, METH_VARARGS,
                   "Call revolve and ask what to do."},
                {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initrevolve(void)
{
  (void) Py_InitModule("revolve", RevolveMethods);
}
