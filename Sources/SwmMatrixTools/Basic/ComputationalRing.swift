//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/10.
//

import SwmCore
#if USE_EIGEN
import SwmEigen
#endif

public protocol ComputationalRing {
    associatedtype ComputationalMatrix: MatrixImpl where ComputationalMatrix.BaseRing == Self
    associatedtype ComputationalSparseMatrix: SparseMatrixImpl where ComputationalSparseMatrix.BaseRing == Self
}

extension Int: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = CSCMatrixImpl<Self>
    #endif
}

extension RationalNumber: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = CSCMatrixImpl<Self>
    #endif
}

extension RealNumber: ComputationalRing {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = CSCMatrixImpl<Self>
}

extension 𝐅₂: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = CSCMatrixImpl<Self>
    #endif
}

extension Polynomial: ComputationalRing where BaseRing: Field {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = CSCMatrixImpl<Self>
}

#if USE_EIGEN
extension EigenSparseMatrixImpl: LUFactorizable where R: EigenSparseMatrixCompatible_LU {}
#endif
