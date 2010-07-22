# -*- coding: utf-8 -*-
import sys
import nuto
import os

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = False

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
A_full_ref = nuto.DoubleFullMatrix(5,5)
A_full_ref.SetValue(0,0,1)
A_full_ref.SetValue(1,0,-2)
A_full_ref.SetValue(3,0,-4)
A_full_ref.SetValue(0,1,-1)
A_full_ref.SetValue(1,1,5)
A_full_ref.SetValue(4,1,8)
A_full_ref.SetValue(2,2,4)
A_full_ref.SetValue(3,2,2)
A_full_ref.SetValue(0,3,-3)
A_full_ref.SetValue(2,3,6)
A_full_ref.SetValue(3,3,7)
A_full_ref.SetValue(2,4,4)
A_full_ref.SetValue(4,4,-5)
if (printResult):
    print "matrix A stored as full matrix (reference)"
    A_full_ref.Info()
    print ""
############################################################################

# sparse matrix with zero based indexing ###################################
A_sparse = nuto.DoubleSparseMatrixCSRGeneral(5,5)
A_sparse.AddEntry(0,0,1)
A_sparse.AddEntry(1,0,-2)
A_sparse.AddEntry(3,0,-4)
A_sparse.AddEntry(0,1,-1)
A_sparse.AddEntry(1,1,5)
A_sparse.AddEntry(4,1,8)
A_sparse.AddEntry(2,2,4)
A_sparse.AddEntry(3,2,2)
A_sparse.AddEntry(0,3,-3)
A_sparse.AddEntry(2,3,6)
A_sparse.AddEntry(3,3,7)
A_sparse.AddEntry(2,4,4)
A_sparse.AddEntry(4,4,-5)
if (printResult):
    print "matrix A, compressed storage, zero based indexing"
    A_sparse.Info()
    print ""

# convert A_sparse into full matrix
A_full = nuto.DoubleFullMatrix(A_sparse)
if (printResult):
    print "matrix A stored as full matrix"
    A_full.Info()
    print ""
A_full -= A_full_ref
if (A_full.Abs().Max()[0] != 0.):
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
B_sparse.AddEntry(0,0,1)
B_sparse.AddEntry(1,0,-2)
B_sparse.AddEntry(3,0,-4)
B_sparse.AddEntry(0,1,-1)
B_sparse.AddEntry(1,1,5)
B_sparse.AddEntry(4,1,8)
B_sparse.AddEntry(2,2,4)
B_sparse.AddEntry(3,2,2)
B_sparse.AddEntry(0,3,-3)
B_sparse.AddEntry(2,3,6)
B_sparse.AddEntry(3,3,7)
B_sparse.AddEntry(2,4,4)
B_sparse.AddEntry(4,4,-5)
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
B_sparse.SetOneBasedIndexing()
B2_sparse.AddEntry(1,1,3)
B2_sparse.AddEntry(0,1,8)
B2_sparse.AddEntry(0,0,2)
B2_sparse.AddEntry(1,3,3)
B2_sparse.AddEntry(1,2,5)
B2_sparse.AddEntry(3,2,1)
B2_sparse.AddEntry(3,3,7)
B2_sparse.AddEntry(3,0,9)
B2_sparse.AddEntry(3,1,4)
if (printResult):
    print "matrix B2, compressed storage, vector of vectors, one based indexing"
    B_sparse.Info()
    print ""
    
if (error):
    sys.exit(-1)
else:
    sys.exit(0)
