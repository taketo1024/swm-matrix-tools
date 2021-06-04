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
    
    // U = [U0, U1]. V = U0^-1 U1.
    // Z = [-V; I] satisfies UZ = O.
    // (PA)(QZ) = (LU)Z = O.
    
    public var kernel: Matrix<m, anySize> {
        let (r, m) = U.size
        let U0 = U.submatrix(colRange: 0 ..< r)
        let U1 = U.submatrix(colRange: r ..< m)
        let V = Matrix.solveUpperTriangular(U0, U1) // U0 * K = U1
        let I = Matrix<anySize, anySize>.identity(size: (m - r, m - r))
        return (-V).stack(I).as(Matrix<m, anySize>.self).permuteRows(by: Q.inverse!)
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
    
    public var cokernelProjector: (Matrix<n, _1>) -> (Matrix<anySize, _1>) {
        let (n, r) = L.size
        let L0 = L.submatrix(rowRange: 0 ..< r)
        
        return { z in
            let Pz = z.permuteRows(by: P)
            let w = Pz.submatrix(rowRange: 0 ..< r)
            let x = Matrix.solveLowerTriangular(L0, w)
            return (Pz - L * x)[r ..< n]
        }
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
        let Pb = b.permuteRows(by: P)
        if let z = Matrix.solveLowerTrapezoidal(L, Pb) {
            let y = Matrix.solveUpperTrapezoidal(U, z)
            let x = y.permuteRows(by: Q.inverse!)
            return x
        } else {
            return nil
        }
    }
    
    // Given z \in Span(Z), solve (QZ)x = z.
    // Put U' = [U0, U1; O, I]. Then
    //
    //   (U'Q^-1)(QZ) = [O; I].
    //
    // so we must have
    //
    //   - U Q^-1 z = 0
    //   - x = (Q^-1 z)[r ..< m].
    //
    public func solveKernel(_ z: ColVectorIF<Impl, m>) -> ColVectorIF<Impl, anySize>? {
        assert(U.size.cols == z.size.rows)
        
        let (r, m) = U.size
        let w = z.permuteRows(by: Q)
        
        if (U * w).isZero {
            return w[r ..< m]
        } else {
            return nil
        }
    }
}
