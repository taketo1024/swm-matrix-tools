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
    associatedtype ComputationalMatrixImpl: MatrixImpl where ComputationalMatrixImpl.BaseRing == Self
    associatedtype ComputationalSparseMatrixImpl: SparseMatrixImpl where ComputationalSparseMatrixImpl.BaseRing == Self
    
    var computationalWeight: Double { get } // used for matrix elimination
}

extension ComputationalRing {
    public typealias ComputationalMatrix<n, m> = MatrixIF<ComputationalMatrixImpl, n, m> where n: SizeType, m: SizeType
    public typealias ComputationalSparseMatrix<n, m> = MatrixIF<ComputationalSparseMatrixImpl, n, m> where n: SizeType, m: SizeType
}

// TODO implement `DefaultDenseMatrixImpl`.

extension Int: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrixImpl = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrixImpl = EigenSparseMatrixImpl<Self>
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
    public typealias ComputationalMatrixImpl = EigenMatrixImpl<Self>
    public typealias ComputationalSparseMatrixImpl = EigenSparseMatrixImpl<Self>
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
    public typealias ComputationalMatrixImpl = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrixImpl = DefaultSparseMatrixImpl<Self>
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double( Swift.abs( max(self, 1/self) ) )
    }
}

extension 𝐅₂: ComputationalRing {
    #if USE_EIGEN
    public typealias ComputationalMatrixImpl = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrixImpl = EigenSparseMatrixImpl<Self>
    #else
    public typealias ComputationalMatrix = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrix = DefaultSparseMatrixImpl<Self>
    #endif
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : 1
    }
}

extension Polynomial: ComputationalRing where BaseRing: Field & ComputationalRing {
    public typealias ComputationalMatrixImpl = DefaultMatrixImpl<Self>
    public typealias ComputationalSparseMatrixImpl = DefaultSparseMatrixImpl<Self>
    
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double(leadExponent + 1) * leadCoeff.computationalWeight
    }
}
