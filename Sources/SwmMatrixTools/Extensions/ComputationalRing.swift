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

// TODO implement `DefaultDenseMatrixImpl`.

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
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
}

extension RationalNumber: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = EigenMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
}

extension RealNumber: ComputationalRing {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
}

extension 𝐅₂: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
}

extension Polynomial: ComputationalRing where BaseRing: Field {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
}
