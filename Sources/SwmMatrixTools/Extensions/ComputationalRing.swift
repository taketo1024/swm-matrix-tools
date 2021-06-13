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
    
    var computationalWeight: Double { get } // used for matrix elimination
}

extension Int: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
    
    @inlinable
    public var computationalWeight: Double {
        Double(Swift.abs(self))
    }
}

extension RationalNumber: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = EigenMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double(max(numerator.abs, denominator))
    }
}

extension RealNumber: ComputationalRing {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double( Swift.abs( max(self, 1/self) ) )
    }
}

extension 𝐅₂: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : 1
    }
}

extension Polynomial: ComputationalRing where BaseRing: Field {
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double(leadExponent + 1)
    }
}
