//
//  LUFactorizer.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/29.
//

import SwmCore

public final class LUFactorizer<M: MatrixImpl & LUFactorizable> where M.BaseRing: EuclideanRing {
    public static func factorize(_ A: M) -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: M, U: M)? {
        let (n, m) = A.size
        let (P1, Q1, L1, U1, S) = partialLU(A) // TODO continue until S is dense enough
        let r1 = L1.size.cols
        
        guard let (P2, Q2, L2, U2) = fullLU(S) else {
            return nil
        }
        
        let r2 = L2.size.cols
        let P = P2.shift(r1) * P1.extend(r2)
        let Q = Q2.shift(r1) * Q1.extend(r2)
        let L = L1.concat(.zero(size: (n - L2.size.rows, r2)).stack(L2))
        let U = U1.stack (.zero(size: (r2, m - U2.size.cols)).concat(L2))
        
        return (P, Q, L, U)
    }
    
    public static func partialLU(_ A: M) -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: M, U: M, S: M) {
        // If
        //
        //   PAQ = [L, B]
        //         [C, D]
        //
        // with L: lower-triangular, then
        //
        //   PAQ = [L, 0] * [I, B']
        //         [C, S]   [0, I ]
        //
        //       = [L] * [I, B'] + [0, 0]
        //         [C]             [0, S] .
        //
        // where
        //
        //    B' = L^-1 B,
        //    S  = D - C L^-1 B.
        //

        let (n, m) = A.size
        let pf = MatrixPivotFinder(A, mode: .colBased)
        pf.run()
        
        let (P, Q) = (pf.rowPermutation, pf.colPermutation)
        let r = pf.pivots.count

        let pAq = A.permute(rowsBy: pf.rowPermutation, colsBy: pf.colPermutation)
        let L = pAq.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)
        let B = pAq.submatrix(rowRange: 0 ..< r, colRange: r ..< m)
        let C = pAq.submatrix(rowRange: r ..< n, colRange: 0 ..< r)
        let D = pAq.submatrix(rowRange: r ..< n, colRange: r ..< m)
        
        let B1 = M.solveLowerTriangular(L, B)!
        let S = D - C * B1
        
        return (
            P: P,
            Q: Q,
            L: L.stack(C),
            U: .identity(size: (r, r)).concat(B1),
            S: S
        )
    }
    
    public static func fullLU(_ A: M) -> (P: Permutation<anySize>, Q: Permutation<anySize>, L: M, U: M)? {
        let e = LUEliminator(data: MatrixEliminationData(A))
        e.run()
        return e.PQLU(M.self)
    }
}

public struct LUFactorization<Impl: MatrixImpl & LUFactorizable, n: SizeType, m: SizeType> {
    public let P: Permutation<n>
    public let Q: Permutation<m>
    public let L: MatrixIF<Impl, n, anySize>
    public let U: MatrixIF<Impl, anySize, m>
    
    public init(P: Permutation<n>, Q: Permutation<m>, L: MatrixIF<Impl, n, anySize>, U: MatrixIF<Impl, anySize, m>) {
        self.P = P
        self.Q = Q
        self.L = L
        self.U = U
    }
    
    public var rank: Int {
        L.size.cols // == U.size.rows
    }
    
    public func solve<k>(_ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k>? {
        if let y = Impl.solveUpperTriangular(U.impl, b.impl),
           let x = Impl.solveLowerTriangular(L.impl, y) {
            return .init(x)
        } else {
            return nil
        }
    }
}

extension Permutation {
    func extend(_ n: Int) -> Permutation<anySize> {
        fatalError("TODO")
    }

    func shift(_ n: Int) -> Permutation<anySize> {
        fatalError("TODO")
    }
}

