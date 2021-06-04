//
//  LUFactorizer.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/29.
//

import SwmCore

public final class LUFactorizer<M: MatrixImpl & LUFactorizable> {
    public typealias Result = (P: Permutation<anySize>, Q: Permutation<anySize>, L: M, U: M)
    
    public static func factorize(_ A: M) -> Result? {
        let (res1, S) = partialLU(A) // TODO continue until S is dense enough
        if S.isZero {
            return res1
        }
        guard let res2 = fullLU(S) else {
            return nil
        }
        return compose(res1, res2)
    }
    
    public static func partialLU(_ A: M) -> (result: Result, S: M) {
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
        
        let B1 = M.solveLowerTriangular(L, B)
        let S = D - C * B1
        
        return (
            result: (
                P: P,
                Q: Q,
                L: L.stack(C),
                U: M.identity(size: (r, r)).concat(B1)
            ),
            S: S
        )
    }
    
    public static func fullLU(_ A: M) -> Result? {
        let e = LUEliminator(data: MatrixEliminationData(A))
        e.run()
        return e.PQLU(M.self)
    }
    
    private static func compose(_ res1: Result, _ res2: Result) -> Result {
        let (P1, Q1, L1, U1) = res1
        let r1 = L1.size.cols

        let (P2s, Q2s, L2s, U2s) = res2
        let r2 = L2s.size.cols
        let P2 = P2s.shifted(r1)
        let Q2 = Q2s.shifted(r1)

        let L2 = M.zero(size: (r1, r2)).stack (L2s)
        let U2 = M.zero(size: (r2, r1)).concat(U2s)
        
        let P = P2 * P1
        let Q = Q2 * Q1
        let L = L1.permuteRows(by: P2).concat(L2)
        let U = U1.permuteCols(by: Q2).stack (U2)
        
        return (P, Q, L, U)
    }
}
