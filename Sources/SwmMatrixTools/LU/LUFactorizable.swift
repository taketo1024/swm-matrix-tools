//
//  LUFactorizable.swift
//  
//
//  Created by Taketo Sano on 2021/06/01.
//

import SwmCore

// MEMO:
// Default implementation for both dense and sparse LU factorization are provided.
// Conforming types can override both or either one to attain better performance.

public protocol LUFactorizable: MatrixImpl {
    typealias RawLUFactorizationResult = (P: Permutation<anySize>, Q: Permutation<anySize>, L: Self, U: Self)
    
    // Find (P, Q, L, U) such that PAQ = LU.
    func LUfactorize() -> RawLUFactorizationResult
    func denseLUfactorize() -> RawLUFactorizationResult
    
    // solve L * x = b.
    static func solveLowerTriangular(_ L: Self, _ b: Self) -> Self
    static func solveLowerTrapezoidal(_ L: Self, _ b: Self) -> Self?

    // solve U * x = b.
    static func solveUpperTriangular(_ U: Self, _ b: Self) -> Self
    static func solveUpperTrapezoidal(_ U: Self, _ b: Self) -> Self
}

extension LUFactorizable {
    public func LUfactorize() -> RawLUFactorizationResult {
        denseLUfactorize()
    }
    
    public func denseLUfactorize() -> RawLUFactorizationResult {
        let f = LUFactorizer(self)
        f.run()
        return f.result
    }
    
    public static func solveLowerTriangular(_ L: Self, _ B: Self) -> Self {
        assert(L.isSquare)
        assert(L.isLowerTriangular)
        assert(L.hasInvertibleDiagonal)
        assert(L.size.rows == B.size.rows)
        
        let (n, k) = B.size
        let entries = Array(0 ..< k).parallelFlatMap { j -> [MatrixEntry<BaseRing>] in
            return self.solveLowerTriangularSingle(L, B, j)
        }
        
        return .init(size: (n, k), entries: entries)
    }
    
    private static func solveLowerTriangularSingle(_ L: Self, _ B: Self, _ j: Int) -> [MatrixEntry<BaseRing>] {
        let n = L.size.rows
        var w = B.colVector(j).serialize()
        
        var x: [MatrixEntry<BaseRing>] = []
        x.reserveCapacity(n)
        
        for i in 0 ..< n {
            let wi = w[i]
            if wi.isZero {
                continue
            }
            
            let l = L.colVector(i)
            let li = l[i, 0]
            let xi = wi * li.inverse!
            
            x.append((i, j, xi))
            
            for (j, _, lj) in l.nonZeroEntries {
                w[j] = w[j] - xi * lj
            }
        }
        
        return x
    }

    public static func solveLowerTrapezoidal(_ L: Self, _ B: Self) -> Self? {
        assert(L.size.rows >= L.size.cols)
        assert(L.size.rows == B.size.rows)
        
        if L.isSquare {
            return solveLowerTriangular(L, B)
        } else {
            let (n, r) = L.size
            let k = B.size.cols
            
            let L0 = L.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)
            let b0 = B.submatrix(rowRange: 0 ..< r, colRange: 0 ..< k)
            let x = solveLowerTriangular(L0, b0)
            
            let L1 = L.submatrix(rowRange: r ..< n, colRange: 0 ..< r)
            let b1 = B.submatrix(rowRange: r ..< n, colRange: 0 ..< k)
            
            if L1 * x == b1 {
                return x
            } else {
                return nil
            }
        }
    }
    
    public static func solveUpperTriangular(_ U: Self, _ B: Self) -> Self {
        assert(U.isSquare)
        assert(U.isUpperTriangular)
        assert(U.hasInvertibleDiagonal)
        assert(U.size.rows == B.size.rows)

        let (n, k) = B.size
        let entries = Array(0 ..< k).parallelFlatMap { j -> [MatrixEntry<BaseRing>] in
            return self.solveUpperTriangularSingle(U, B, j)
        }
        
        return .init(size: (n, k), entries: entries)
    }
    
    private static func solveUpperTriangularSingle(_ U: Self, _ B: Self, _ j: Int) -> [MatrixEntry<BaseRing>] {
        let n = U.size.rows
        var w = B.colVector(j).serialize()
        
        var x: [MatrixEntry<BaseRing>] = []
        x.reserveCapacity(n)
        
        for i in (0 ..< n).reversed() {
            let wi = w[i]
            if wi.isZero {
                continue
            }

            let u = U.colVector(i)
            let ui = u[i, 0]
            let xi = wi * ui.inverse!

            x.append((i, j, xi))
            
            for (j, _, uj) in u.nonZeroEntries {
                w[j] = w[j] - xi * uj
            }
        }
        
        return x
    }

    public static func solveUpperTrapezoidal(_ U: Self, _ b: Self) -> Self {
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

public protocol SparseLUFactorizable: SparseMatrixImpl, LUFactorizable {
    func sparseLUfactorize() -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: Self, U: Self)
}

extension LUFactorizable where Self: SparseMatrixImpl {
    public func LUfactorize() -> RawLUFactorizationResult {
        sparseLUfactorize()
    }
    
    public func sparseLUfactorize() -> RawLUFactorizationResult {
        let f = SparseLUFactorizer(self)
        f.run()
        return f.result
    }
}
