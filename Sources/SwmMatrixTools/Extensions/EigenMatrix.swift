//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/12.
//

#if USE_EIGEN
import SwmEigen

extension EigenMatrixImpl: LUFactorizable where R: ComputationalRing & EigenMatrixCompatible_LU {}
extension EigenSparseMatrixImpl: LUFactorizable where R: ComputationalRing & EigenSparseMatrixCompatible_LU {}

#endif

