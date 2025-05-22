import numpy as np
cimport numpy as np

from cpython.ref cimport Py_INCREF

from cython.operator cimport dereference as deref
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np
from libc.math cimport sqrt

np.import_array()

# Context:
#
#   In some Cython implementations, one may need resizable buffers
#   which they would like to return as np.ndarray to Python callers.
#
#   In term of technical choice, C++ std::vector comes as handy
#   datastructures in Cython for interacting efficiently with resizable
#   buffers.
#
#   Yet to my knowledge, there isn't a way to return those vectors to
#   python callers within Cython, especially as np.ndarray.
#
#   To return numpy arrays wrapping those resizable buffers,
#   we can use PyArray_SimpleNewFromData by providing a
#   pointer to the start of the buffer with some metadata.
#
#   Though it is possible to access the buffer first element address
#   with std::vector::data, the buffer themselves can't be stolen:
#   the buffers' lifetime is tight to their std::vectors and are
#   deallocated when their std::vectors are.
#
# Solution proposal:
#
#   To solve this, we propose introducing a StdVectorSentinel, a
#   sentinel to be used as np.ndarrays' base objects, allowing
#   performing a proper coercion between std::vector and np.ndarrays.
#
#   This way, std::vectors can dynamically be allocated in Cython
#   and returned as numpy arrays without inconsistency nor memory leaks.
#
# Question and further work:
#
#   1. Is there already a way to coerce std::vector in numpy array efficiently
#   within Cython?
#   2. If not, is the proposed fixture this a good technical solution?
#   3. If so, could it further be improved to support type covariances, allowing
#   using it easily dtype-invariantly?

####

# The following defines portable dtypes that are to be used for vectors.
cdef enum:
    DTYPECODE = np.NPY_FLOAT64
    ITYPECODE = np.NPY_INTP

ctypedef np.float64_t DTYPE_t
ctypedef np.intp_t ITYPE_t

ctypedef fused DITYPE_t:
    ITYPE_t
    DTYPE_t

ITYPE = np.intp
DTYPE = np.float64

## std::vector to np.ndarray coercion
#
#   As type covariance is not supported for C++ containers via Cython,
#   we need to redefine fused types for vectors.
ctypedef fused vector_DITYPE_t:
    vector[ITYPE_t]
    vector[DTYPE_t]


ctypedef fused vector_vector_DITYPE_t:
    vector[vector[ITYPE_t]]
    vector[vector[DTYPE_t]]


cdef class StdVectorSentinel:
    """Wraps a reference to a vector which will be deallocated with this object.
    When created, the StdVectorSentinel swaps the reference of its internal
    vectors with the provided one (vec_ptr), thus making the StdVectorSentinel
    manage the provided one's lifetime.
    """
    pass


# We necessarily need to define two extension types extending StdVectorSentinel
# because we need to provide the dtype of the vector but can't use numeric fused types.
cdef class StdVectorSentinelDTYPE(StdVectorSentinel):
    cdef vector[DTYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[DTYPE_t] * vec_ptr):
        # This initializes the object directly without calling __init__
        cdef StdVectorSentinelDTYPE sentinel = StdVectorSentinelDTYPE.__new__(StdVectorSentinelDTYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef class StdVectorSentinelITYPE(StdVectorSentinel):
    cdef vector[ITYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[ITYPE_t] * vec_ptr):
        # This initializes the object directly without calling __init__
        cdef StdVectorSentinelITYPE sentinel = StdVectorSentinelITYPE.__new__(StdVectorSentinelITYPE)
        sentinel.vec.swap(deref(vec_ptr))
        return sentinel


cdef np.ndarray vector_to_nd_array(vector_DITYPE_t * vect_ptr):
    """Create a numpy ndarray given a C++ vector.
    The numpy array buffer is the one of the C++ vector.
    A StdVectorSentinel is registered as the base object for the numpy array,
    freeing the C++ vector it encapsulates when the numpy array is freed.
    """
    typenum = DTYPECODE if vector_DITYPE_t is vector[DTYPE_t] else ITYPECODE
    cdef:
        np.npy_intp size = deref(vect_ptr).size()
        np.ndarray arr = np.PyArray_SimpleNewFromData(1, &size, typenum,
                                                      deref(vect_ptr).data())
        StdVectorSentinel sentinel

    if vector_DITYPE_t is vector[DTYPE_t]:
        sentinel = StdVectorSentinelDTYPE.create_for(vect_ptr)
    else:
        sentinel = StdVectorSentinelITYPE.create_for(vect_ptr)

    # Makes the numpy array responsible of the life-cycle of its buffer.
    # A reference to the StdVectorSentinel will be stolen by the call bellow,
    # so we increase its reference counter.
    # See: https://docs.python.org/3/c-api/intro.html#reference-count-details
    Py_INCREF(sentinel)
    np.PyArray_SetBaseObject(arr, sentinel)
    return arr

# Module callable front-end for interactive use and proof of concept.
cpdef np.ndarray get_np_array_std_vector_backed(ITYPE_t n_elements):

    cdef vector[ITYPE_t] * vector_of_ints = new vector[ITYPE_t](n_elements, 1337)

    return vector_to_nd_array(vector_of_ints)