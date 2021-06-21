import os
from rpy2.robjects import r, numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
numpy2ri.activate()
pandas2ri.activate()
importr('gifti')

neurobase = importr('neurobase')

def readnii(file):
	return neurobase.readnii(file)
def writenii(data, file):
	import nibabel as nib
	import numpy as np
	if type(data) is np.ndarray:
		new_image = nib.Nifti1Image(data, affine=np.eye(4))
		return new_image.to_filename(file)
	else:
		return neurobase.writenii(data, file)


def convert_binary_datatype(datatype = None):
	
	if datatype is None:
		datatype = r('c("NIFTI_TYPE_UINT8", "NIFTI_TYPE_INT32", "NIFTI_TYPE_UINT32",
    "NIFTI_TYPE_FLOAT32")')
	r.assign('datatype', datatype)
	return r('convert_binary_datatype(datatype=datatype)')

def convert_endian(endian = r('NA')):
	r.assign('endian', endian)
	return r('convert_endian(endian=endian)')

def convert_intent(intent = r('NA')):
	r.assign('intent', intent)
	return r('convert_intent(intent=intent)')

def create_data_matrix(data = r('NA'), dims = r('NA'), ordering = None):
	r.assign('data', data)
	r.assign('dims', dims)
	if ordering is None:
		ordering = r('c("RowMajorOrder", "ColumnMajorOrder")')
	r.assign('ordering', ordering)
	return r('create_data_matrix(data=data, dims=dims, ordering=ordering)')

def data_array_attributes(darray = r('NA')):
	r.assign('darray', darray)
	return r('data_array_attributes(darray=darray)')

def data_decoder(values = r('NA'), encoding = None, datatype = r('NULL'), endian = None, ext_filename = r('NULL'), n = None):
	r.assign('values', values)
	r.assign('datatype', datatype)
	r.assign('ext_filename', ext_filename)
	if encoding is None:
		encoding = r('c("ASCII", "Base64Binary", "GZipBase64Binary", "ExternalFileBinary")')
	r.assign('encoding', encoding)
	if endian is None:
		endian = r('c("little", "big", "LittleEndian", "BigEndian")')
	r.assign('endian', endian)
	if n is None:
		n = r('defi')
	r.assign('n', n)
	return r('data_decoder(values=values, encoding=encoding, datatype=datatype, endian=endian, ext_filename=ext_filename, n=n)')

def data_encoder(values = r('NA'), encoding = None, datatype = r('NULL'), endian = None):
	r.assign('values', values)
	r.assign('datatype', datatype)
	if encoding is None:
		encoding = r('c("ASCII", "Base64Binary", "GZipBase64Binary")')
	r.assign('encoding', encoding)
	if endian is None:
		endian = r('c("little", "big", "LittleEndian", "BigEndian")')
	r.assign('endian', endian)
	return r('data_encoder(values=values, encoding=encoding, datatype=datatype, endian=endian)')

def decompress_gii(file = r('NA')):
	r.assign('file', file)
	return r('decompress_gii(file=file)')

def download_gifti_data(outdir = None, overwrite = False):
	r.assign('overwrite', overwrite)
	if outdir is None:
		outdir = r('system.file(package = "gifti")')
	r.assign('outdir', outdir)
	return r('download_gifti_data(outdir=outdir, overwrite=overwrite)')

def gifti_list(file = r('NA')):
	r.assign('file', file)
	return r('gifti_list(file=file)')

def gifti_map_value(pointset = r('NA'), triangle = r('NA'), values = r('NA'), indices = None, add_one = True):
	r.assign('pointset', pointset)
	r.assign('triangle', triangle)
	r.assign('values', values)
	r.assign('add_one', add_one)
	if indices is None:
		indices = r('seq(nrow(pointset))')
	r.assign('indices', indices)
	return r('gifti_map_value(pointset=pointset, triangle=triangle, values=values, indices=indices, add_one=add_one)')

def have_gifti_test_data(outdir = None):
	
	if outdir is None:
		outdir = r('system.file(package = "gifti"')
	r.assign('outdir', outdir)
	return r('have_gifti_test_data(outdir=outdir)')

def is_gifti(x = None):
	
	if x is None:
		x = r(')

is_gifti(')
	r.assign('x', x)
	return r('is.gifti(x=x)')

def readgii(file = None):
	
	if file is None:
		file = r(')

readGIfTI(')
	r.assign('file', file)
	return r('readgii(file=file)')

def surf_triangles(file = r('NA')):
	r.assign('file', file)
	return r('surf_triangles(file=file)')

def writegii(gii = r('NA'), out_file = r('NA'), use_parsed_transformations = None):
	r.assign('gii', gii)
	r.assign('out_file', out_file)
	if use_parsed_transformations is None:
		use_parsed_transformations = r('FALSE)

writeGIfTI(gii, out_file')
	r.assign('use_parsed_transformations', use_parsed_transformations)
	return r('writegii(gii=gii, out_file=out_file, use_parsed_transformations=use_parsed_transformations)')
