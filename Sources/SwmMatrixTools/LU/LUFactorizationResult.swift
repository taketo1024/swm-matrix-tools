//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/02.
//

import SwmCore

public struct LUFactorizationResult<Impl: MatrixImpl & LUFactorizable, n: SizeType, m: SizeType> {
    public typealias Matrix<n, m> = MatrixIF<Impl, n, m> where n: SizeType, m: SizeType
    public typealias MatrixL = Matrix<n, anySize>
    public typealias MatrixU = Matrix<anySize, m>
    public typealias PermutationP  = Permutation<n>
    public typealias PermutationQ = Permutation<m>

    public let P: PermutationP
    public let Q: PermutationQ
    public let L: MatrixL
    public let U: MatrixU

    public init(P: Permutation<n>, Q: Permutation<m>, L: MatrixL, U: MatrixU) {
        self.P = P
        self.Q = Q
        self.L = L
        self.U = U
    }
    
    public var PQLU: (PermutationP, PermutationQ, MatrixL, MatrixU) {
        (P, Q, L, U)
    }
    
    public var rank: Int {
        L.size.cols
    }
    
    public var nullity: Int {
        U.size.cols - rank
    }
    
    public var kernel: Matrix<m, anySize> {
        let (r, m) = U.size
        let U0 = U.submatrix(colRange: 0 ..< r)
        let U1 = U.submatrix(colRange: r ..< m)
        let K = Matrix.solveUpperTriangular(U0, U1)! // U0 * K = U1
        let I = Matrix<anySize, anySize>.identity(size: (m - r, m - r))
        return (-K).stack(I).as(Matrix<m, anySize>.self).permuteRows(by: Q.inverse!)
    }
    
    // Im(A) = Im(P^-1 L U Q^-1) = Im(P^-1 L).
    public var image: Matrix<n, anySize> {
        L.permuteRows(by: P.inverse!)
    }
    
    // Coker(f) := W / Im(f) = Im(P^-1 [O, I_{n-r}])
    public var cokernel: Matrix<n, anySize> {
        let (n, r) = L.size
        let I = Matrix<n, anySize>.colUnits(size: (n, n - r), indices: (r ..< n))
        return I.permuteRows(by: P.inverse!)
    }
    
    public func solve<k>(_ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k>? {
        // Solve Ax = b
        //
        //  <=> PAx = (PAQ)(Q^-1 x)
        //          = LU(Q^-1 x)
        //          = Pb.
        //
        //  <=> Lz = Pb, Uy = z, x = Qy.
        
        assert(L.size.rows == b.size.rows)
        let Pb = b.permuteRows(by: P).impl
        if let z = Impl.solveLowerTriangular(L.impl, Pb),
           let y = Impl.solveUpperTriangular(U.impl, z) {
            //  Recall Qy = (Q1 ... Qn)y
            let x = y.permuteRows(by: Q.inverse!.asAnySize)
            return .init(x)
        } else {
            return nil
        }
    }
}
