//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/02.
//

import SwmCore

extension MatrixImpl {
    public var isUpperTriangular: Bool {
        nonZeroEntries.allSatisfy{ (i, j, a) in i <= j }
    }
    
    public var isLowerTriangular: Bool {
        nonZeroEntries.allSatisfy{ (i, j, a) in i >= j }
    }
    
    public var hasInvertibleDiagonal: Bool {
        nonZeroEntries.allSatisfy{ (i, j, a) in
            i != j || a.isInvertible
        }
    }
}
