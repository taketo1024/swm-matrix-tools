//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/09/30.
//

import SwmCore

extension DefaultSparseMatrixImpl where BaseRing: Field {
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
}

extension DefaultSparseMatrixImpl: LUFactorizable where BaseRing: Field & ComputationalRing {
    public func LUfactorize() -> RawLUFactorizationResult {
        let f = SparseLUFactorizer(self)
        f.run()
        return f.result
    }
}
