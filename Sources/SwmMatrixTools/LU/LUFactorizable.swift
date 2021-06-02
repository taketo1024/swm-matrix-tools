//
//  LUFactorizable.swift
//  
//
//  Created by Taketo Sano on 2021/06/01.
//

import SwmCore

public protocol LUFactorizable: MatrixImpl {
    // solve L * x = b.
    static func solveLowerTriangular(_ L: Self, _ b: Self) -> Self

    // solve U * x = b.
    static func solveUpperTriangular(_ U: Self, _ b: Self) -> Self
}

extension LUFactorizable {
    static func solveLowerTrapezoidal(_ L: Self, _ b: Self) -> Self? {
        assert(L.size.rows >= L.size.cols)
        assert(L.size.rows == b.size.rows)
        
        if L.isSquare {
            return solveLowerTriangular(L, b)
        } else {
            let (n, r) = L.size
            let k = b.size.cols
            
            let L0 = L.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)
            let b0 = b.submatrix(rowRange: 0 ..< r, colRange: 0 ..< k)
            let x = solveLowerTriangular(L0, b0)
            
            let L1 = L.submatrix(rowRange: r ..< n, colRange: 0 ..< r)
            let b1 = b.submatrix(rowRange: r ..< n, colRange: 0 ..< k)
            
            if L1 * x == b1 {
                return x
            } else {
                return nil
            }
        }
    }
    
    static func solveUpperTrapezoidal(_ U: Self, _ b: Self) -> Self {
        assert(U.size.rows <= U.size.cols)
        assert(U.size.rows == b.size.rows)

        if U.isSquare {
            return solveUpperTriangular(U, b)
        } else {
            let (r, m) = U.size
            let k = b.size.cols
            
            let U0 = U.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)
            let x = solveUpperTriangular(U0, b)
            
            return x.stack(.zero(size: (m - r, k)))
        }
    }
}
