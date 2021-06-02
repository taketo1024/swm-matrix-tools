//
//  LUFactorizable.swift
//  
//
//  Created by Taketo Sano on 2021/06/01.
//

import SwmCore

public protocol LUFactorizable: MatrixImpl {
    func factorize() -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: Self, U: Self)
    
    // solve L * x = b.
    static func solveLowerTriangular(_ L: Self, _ b: Self) -> Self?
    
    // solve U * x = b.
    static func solveUpperTriangular(_ U: Self, _ b: Self) -> Self?
}

extension LUFactorizable {
    public static func solveLowerTriangular(_ L: Self, _ B: Self) -> Self? {
        assert(L.isLowerTriangular && L.hasInvertibleDiagonal)
        
        let (m, k) = (L.size.cols, B.size.cols)
        let Lcols = L.splitIntoCols()
        let Bcols = B.splitIntoCols()
        let Xcols = Bcols.parallelMap { b in
            self.solveLowerTriangularSingle(Lcols, b)
        }
        
        if Xcols.contains(where: { $0 == nil }) {
            return nil
        } else {
            return .init(size: (m, k), cols: Xcols.map{ $0! } )
        }
    }
    
    public static func solveLowerTriangularSingle(_ Lcols: [Self], _ b: Self) -> Self? {
        typealias R = BaseRing
        
        let n = b.size.rows
        let m = Lcols.count
        let r = min(n, m)
        
        var x = Self.zero(size: (m, 1))
        var w = b
        
        for i in 0 ..< r {
            let wi = w[i, 0]
            if wi.isZero {
                continue
            }
            
            let l = Lcols[i]
            let e = l[i, 0]
            let xi = wi * e.inverse!
            
            x[i, 0] = xi
            w = w - xi * l
        }
        if w.isZero {
            return x
        } else {
            return nil
        }
    }

    public static func solveUpperTriangular(_ U: Self, _ B: Self) -> Self? {
        assert(U.isUpperTriangular && U.hasInvertibleDiagonal)
        
        let (n, m) = U.size
        let r = min(n, m)
        let k = B.size.cols
        
        if !B.submatrix(rowRange: r ..< n, colRange: 0 ..< k).isZero {
            return nil
        }
        
        let Ucols = U.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r).splitIntoCols()
        let Bcols = B.submatrix(rowRange: 0 ..< r, colRange: 0 ..< k).splitIntoCols()
        let Xcols = Bcols.parallelMap { b -> Self in
            let x = self.solveUpperTriangularSingle(Ucols, b)
            return x.stack(.zero(size: (n - r, 1)))
        }
        return .init(size: (m, k), cols: Xcols)
    }
    
    private static func solveUpperTriangularSingle(_ Ucols: [Self], _ b: Self) -> Self {
        let r = b.size.rows
        var x = Self.zero(size: (r, 1))
        var w = b
        
        for i in (0 ..< r).reversed() {
            let wi = w[i, 0]
            if wi.isZero {
                continue
            }
            
            let u = Ucols[i]
            let e = u[i, 0]
            let xi = wi * e.inverse!
            
            x[i, 0] = xi
            w = w - xi * u
        }
        return x
    }
    
    private func splitIntoCols() -> [Self] {
        let cols = nonZeroEntries.group { $0.col }
        return (0 ..< size.cols).map { j in
            let col = cols[j] ?? []
            return Self(
                    size: (size.rows, 1),
                    entries: col.map{ (i, _, a) in MatrixEntry(i, 0, a) }
                )
            }
    }
    
    private init(size: MatrixSize, cols: [Self]) {
        self.init(
            size: size,
            entries: cols.enumerated().flatMap { (j, x) in
                x.nonZeroEntries.map { (i, _, a) in
                    (i, j, a)
                }
            }
        )
    }
}

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

extension DefaultMatrixImpl: LUFactorizable where BaseRing: Field {
    public func factorize() -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: DefaultMatrixImpl<R>, U: DefaultMatrixImpl<R>) {
        fatalError("")
    }
}
