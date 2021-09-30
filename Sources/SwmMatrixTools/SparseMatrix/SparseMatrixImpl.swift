//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/09/30.
//

import SwmCore

public protocol SparseMatrixImpl: MatrixImpl {
    var numberOfNonZeros: Int { get }
    var density: Double { get }
}

extension SparseMatrixImpl {
    public var isZero: Bool {
        numberOfNonZeros == 0
    }
    
    public var density: Double {
        let N = numberOfNonZeros
        return N > 0 ? Double(N) / Double(size.rows * size.cols) : 0
    }
}
