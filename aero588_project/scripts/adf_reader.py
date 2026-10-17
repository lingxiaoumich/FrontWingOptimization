"""Read an ADF CGNS tree through the official CGNS ADF library.

The library must use the CGNS 4.x ADF API and 64-bit cgsize_t. Opened databases
are read-only. Pass a suitable shared-library path; no file data is modified.
"""
import ctypes as ct
import numpy as np

class ADFReader:
    def __init__(self, library):
        self.lib=ct.CDLL(str(library))

    def call(self, name, *args):
        error=ct.c_int()
        getattr(self.lib,'ADF_'+name)(*args,ct.byref(error))
        if error.value != -1:
            message=ct.create_string_buffer(100)
            self.lib.ADF_Error_Message(error,message)
            raise RuntimeError(name+': '+message.value.decode())

    def text(self, node, field):
        buffer=ct.create_string_buffer(100)
        self.call('Get_'+field,ct.c_double(node),buffer)
        return buffer.value.decode().strip()

    def children(self, node):
        n=ct.c_int()
        self.call('Number_of_Children',ct.c_double(node),ct.byref(n))
        ids=(ct.c_double*n.value)()
        count=ct.c_int()
        if n.value:
            self.call('Children_IDs',ct.c_double(node),ct.c_int(1),n,ct.byref(count),ids)
        return list(ids)

    def read(self, node):
        dtype=self.text(node,'Data_Type')
        if dtype=='MT':
            return None
        ndims=ct.c_int()
        self.call('Get_Number_of_Dimensions',ct.c_double(node),ct.byref(ndims))
        dims=(ct.c_longlong*12)()
        self.call('Get_Dimension_Values',ct.c_double(node),dims)
        shape=tuple(dims[:ndims.value])
        types={'R4':'float32','R8':'float64','I4':'int32','I8':'int64','C1':'S1'}
        array=np.empty(shape,dtype=types[dtype],order='F')
        self.call('Read_All_Data',ct.c_double(node),dtype.encode(),ct.c_void_p(array.ctypes.data))
        return array

    def load(self, path):
        root=ct.c_double()
        self.call('Database_Open',str(path).encode(),b'READ_ONLY',b'NATIVE',ct.byref(root))
        def recurse(node):
            return {'name':self.text(node,'Name'),'label':self.text(node,'Label'),
                    'data':self.read(node),'children':[recurse(c) for c in self.children(node)]}
        try:
            return recurse(root.value)
        finally:
            self.call('Database_Close',root)
