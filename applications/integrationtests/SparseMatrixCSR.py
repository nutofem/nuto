# -*- coding: utf-8 -*-
import sys
import nuto
import os
import numpy as np

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = True

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% start the real test file                                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create a full matrix as reference ########################################
A_full_ref = np.zeros((5,5))
A_full_ref[0,0] = 1.
A_full_ref[1,0] = -2.
A_full_ref[3,0] = -4.
A_full_ref[0,1] = -1.
A_full_ref[1,1] = 5.
A_full_ref[4,1] = 8.
A_full_ref[2,2] = 4.
A_full_ref[3,2] = 2.
A_full_ref[0,3] = -3.
A_full_ref[2,3] = 6.
A_full_ref[3,3] = 7.
A_full_ref[2,4] = 4.
A_full_ref[4,4] = -5.
if (printResult):
    print "matrix A stored as full matrix (reference)"
    print A_full_ref
    print ""
############################################################################

# sparse matrix with zero based indexing ###################################
A_sparse = nuto.DoubleSparseMatrixCSRGeneral(5,5)
A_sparse.AddValue(0,0,1)
A_sparse.AddValue(1,0,-2)
A_sparse.AddValue(3,0,-4)
A_sparse.AddValue(0,1,-1)
A_sparse.AddValue(1,1,5)
A_sparse.AddValue(4,1,8)
A_sparse.AddValue(2,2,4)
A_sparse.AddValue(3,2,2)
A_sparse.AddValue(0,3,-3)
A_sparse.AddValue(2,3,6)
A_sparse.AddValue(3,3,7)
A_sparse.AddValue(2,4,4)
A_sparse.AddValue(4,4,-5)
if (printResult):
    print "matrix A, compressed storage, zero based indexing"
    A_sparse.Info()
    print ""

# convert A_sparse into full matrix
A_full = A_sparse.ConvertToFullMatrix()
if (printResult):
    print "matrix A stored as full matrix"
    print A_full
    print ""
A_full -= A_full_ref
if (np.max(np.abs(A_full)) != 0.):
    print '[' + system,sys.argv[0] + '] : error converting sparse matrix with zero indexing to full matrix.'
    error = True;

#if create_reference:
#    A_sparse.save(path_xml_sources+"sparseMatrixCSR_result_1_"+system+".xml",nuto.XML)
#else:
#    A_sparse_ref = nuto.DoubleSparseMatrixCSRGeneral(0,0)
#    A_sparse_ref.restore(path_xml_sources+"sparseMatrixCSR_result_1_"+system+".xml",nuto.XML)

# switch to one based indexing
A_sparse.SetOneBasedIndexing()
if (printResult):
    print "matrix a, one based indexing"
    A_sparse.Info()
    print ""
#if create_reference:
#    A_sparse.save(path_xml_sources+"sparseMatrixCSR_result_2_"+system+".xml",nuto.XML)
#else:
#    A_sparse_ref = nuto.DoubleSparseMatrixCSRGeneral(0,0)
#    A_sparse_ref.restore(path_xml_sources+"sparseMatrixCSR_result_2_"+system+".xml",nuto.XML)

############################################################################
# sparse matrix with one based indexing ####################################
# (interface still uses zero based indexing)
B_sparse = nuto.DoubleSparseMatrixCSRGeneral(5,5,13)
B_sparse.SetOneBasedIndexing()
B_sparse.AddValue(0,0,1)
B_sparse.AddValue(1,0,-2)
B_sparse.AddValue(3,0,-4)
B_sparse.AddValue(0,1,-1)
B_sparse.AddValue(1,1,5)
B_sparse.AddValue(4,1,8)
B_sparse.AddValue(2,2,4)
B_sparse.AddValue(3,2,2)
B_sparse.AddValue(0,3,-3)
B_sparse.AddValue(2,3,6)
B_sparse.AddValue(3,3,7)
B_sparse.AddValue(2,4,4)
B_sparse.AddValue(4,4,-5)
if (printResult):
    print "matrix B, compressed storage, one based indexing"
    B_sparse.Info()
    print ""
#if create_reference:
#    B_sparse.save(path_xml_sources+"sparseMatrixCSR_result_3_"+system+".xml",nuto.XML)
#else:
#    B_sparse_ref = nuto.DoubleSparseMatrixCSRGeneral(0,0)
#    B_sparse_ref.restore(path_xml_sources+"sparseMatrixCSR_result_3_"+system+".xml",nuto.XML)

# switch to zero based indexing
B_sparse.SetZeroBasedIndexing()
if (printResult):
    print "matrix B, compressed storage, zero based indexing"
    B_sparse.Info()
    print ""
#if create_reference:
#    B_sparse.save(path_xml_sources+"sparseMatrixCSR_result_4_"+system+".xml",nuto.XML)
#else:
#    B_sparse_ref = nuto.DoubleSparseMatrixCSRGeneral(0,0)
#    B_sparse_ref.restore(path_xml_sources+"sparseMatrixCSR_result_4_"+system+".xml",nuto.XML)

############################################################################
# sparse matrix vector of vector with one based indexing ###################
# (interface still uses zero based indexing)
B2_sparse = nuto.DoubleSparseMatrixCSRVector2General(4,4)
B2_sparse.AddValue(1,1,3)
B2_sparse.AddValue(0,1,8)
B2_sparse.AddValue(0,0,2)
B2_sparse.AddValue(1,3,3)
B2_sparse.AddValue(1,2,5)
B2_sparse.AddValue(3,2,1)
B2_sparse.AddValue(3,3,7)
B2_sparse.AddValue(3,0,9)
B2_sparse.AddValue(3,1,4)
if (printResult):
    print "matrix B2, compressed storage, vector of vectors, one based indexing"
    B_sparse.Info()
    print "matrix B2, converted to full matrix"
    print B2_sparse.ConvertToFullMatrix()
    print ""

C2_sparse = nuto.DoubleSparseMatrixCSRVector2General(4,1)
C2_sparse.AddValue(0,0,4)
C2_sparse.AddValue(1,0,2)
C2_sparse.AddValue(3,0,3)
if (printResult):
    print "matrix C2, converted to full matrix"
    print C2_sparse.ConvertToFullMatrix()
    print ""

D2_sparse = nuto.DoubleSparseMatrixCSRVector2General(1,5)
D2_sparse.AddValue(0,1,4)
D2_sparse.AddValue(0,4,1)
if (printResult):
    print "matrix D2, converted to full matrix"
    print D2_sparse.ConvertToFullMatrix()
    print ""

B2_sparse.ConcatenateColumns(C2_sparse)
B2_sparse.ConcatenateRows(D2_sparse)
B2_Full = B2_sparse.ConvertToFullMatrix()
if (printResult):
    print "matrix B2, appended Columns of C2 and then append rows of D2"
    print B2_Full
    print ""
B2_FullRef = np.zeros((5,5))
B2_FullRef[1,1] = 3.
B2_FullRef[0,1] = 8.
B2_FullRef[0,0] = 2.
B2_FullRef[1,3] = 3.
B2_FullRef[1,2] = 5.
B2_FullRef[3,2] = 1.
B2_FullRef[3,3] = 7.
B2_FullRef[3,0] = 9.
B2_FullRef[3,1] = 4.
B2_FullRef[0,4] = 4.
B2_FullRef[1,4] = 2.
B2_FullRef[3,4] = 3.
B2_FullRef[4,1] = 4.
B2_FullRef[4,4] = 1.
if (printResult):
    print "matrix B2_ref"
    print B2_FullRef
    print ""

if (np.max(np.abs(B2_Full-B2_FullRef)) > 1e-8):
    print '[' + system,sys.argv[0] + '] : concatenation of B2 is not correct.'
    error = True
    
if (error):
    sys.exit(-1)
else:
    sys.exit(0)
