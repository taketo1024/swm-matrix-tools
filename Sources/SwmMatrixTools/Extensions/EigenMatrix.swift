//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/12.
//

#if USE_EIGEN
import SwmEigen

extension EigenSparseMatrixImpl: LUFactorizable, SparseLUFactorizable where R: EigenSparseMatrixCompatible_LU {}

#endif

