//
//  MatrixEliminationResult.swift
//  Sample
//
//  Created by Taketo Sano on 2018/04/26.
//

import SwmCore

public struct MatrixEliminationResult<Impl: MatrixImpl, n: SizeType, m: SizeType> {
    public typealias R = Impl.BaseRing
    
    public let form: MatrixEliminationForm
    public let size: MatrixSize
    public let entries: AnySequence<MatrixEntry<R>>
    public let headEntries: AnySequence<MatrixEntry<R>>
    public let rowOps: [RowElementaryOperation<R>]
    public let colOps: [ColElementaryOperation<R>]
    
    public init(form: MatrixEliminationForm, size: MatrixSize, entries: AnySequence<MatrixEntry<R>>, headEntries: AnySequence<MatrixEntry<R>>, rowOps: [RowElementaryOperation<R>], colOps: [ColElementaryOperation<R>]) {
        self.form = form
        self.size = size
        self.entries = entries
        self.headEntries = headEntries
        self.rowOps = rowOps
        self.colOps = colOps
    }
    
    public var transposed: MatrixEliminationResult<Impl, m, n> {
        func transpose<S: Sequence>(_ entries: S) -> AnySequence<MatrixEntry<R>> where S.Element == MatrixEntry<R> {
            AnySequence(entries.lazy.map{ (i, j, a) in (j, i, a) })
        }
        return .init(
            form: form.transposed,
            size: (size.cols, size.rows),
            entries: transpose(entries),
            headEntries: transpose(headEntries),
            rowOps: colOps.map{ $0.transposed },
            colOps: rowOps.map{ $0.transposed }
        )
    }
    
    public var result: MatrixIF<Impl, n, m> {
        .init(size: size, entries: entries)
    }
    
    public var isSquare: Bool {
        size.rows == size.cols
    }
    
    public var isDiagonal: Bool {
        entries.allSatisfy { (i, j, _) in i == j }
    }
    
    public var isIdentity: Bool {
        isSquare && entries.allSatisfy { (i, j, a) in i == j && a.isIdentity }
    }
    
    public var rank: Int {
        assert(form != .none)
        return headEntries.count
    }
    
    public var nullity: Int {
        size.cols - rank
    }
    
    public var left: MatrixIF<Impl, n, n> {
        .identity(size: (size.rows, size.rows))
            .appliedRowOperations(rowOps)
    }
    
    public var leftInverse: MatrixIF<Impl, n, n> {
        .identity(size: (size.rows, size.rows))
            .appliedRowOperations(rowOpsInverse)
    }
    
    public var right: MatrixIF<Impl, m, m> {
        .identity(size: (size.cols, size.cols))
            .appliedColOperations(colOps)
    }
    
    public var rightInverse: MatrixIF<Impl, m, m> {
        .identity(size: (size.cols, size.cols))
            .appliedColOperations(colOpsInverse)
    }
    
    public var rowOpsInverse: [RowElementaryOperation<R>] {
        rowOps.reversed().map{ $0.inverse }
    }

    public var colOpsInverse: [ColElementaryOperation<R>] {
        colOps.reversed().map{ $0.inverse }
    }
}

extension MatrixEliminationResult where n == m {
    public var determinant: R {
        assert(isSquare)
        if rank == size.rows {
            return rowOps.multiply { $0.determinant }.inverse!
                * colOps.multiply { $0.determinant }.inverse!
                * headEntries.map{ $0.value }.multiply()
        } else {
            return .zero
        }
    }
    
    public var inverse: MatrixIF<Impl, n, n>? {
        assert(isSquare)
        return (isIdentity) ? right * left : nil
    }
}

extension MatrixEliminationResult where R: EuclideanRing {
    // eliminate again
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
        let e = result.eliminate(form: form)
        return .init(
            form: form,
            size: size,
            entries: e.entries,
            headEntries: e.headEntries,
            rowOps: rowOps + e.rowOps,
            colOps: colOps + e.colOps
        )
    }
}

// MARK: Associated matrices

extension MatrixEliminationResult {
    
    // Returns the matrix B consisting of the basis vectors of Im(A).
    // With
    //
    //     P * A * Q = [ L  | O   ]
    //                 [ L' | O_k ]
    //
    // Im(A) is spanned by [L; *].
    //
    //  => B = P^-1 * [L; L'].
    //
    
    public var imageMatrix: MatrixIF<Impl, n, anySize> {
        let (r, n) = (rank, size.rows)
        if r == 0 {
            return .zero(size: (n, 0))
        }
        
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return result
                .submatrix(colRange: 0 ..< r)
                .appliedRowOperations(rowOpsInverse)

        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // Returns a matrix Z consisting of the basis vectors of Ker(A).
    // With k = m - rank,
    //
    //     P * A * Q = [ L  | O   ]
    //                 [ L' | O_k ]
    //
    //  => Z = Q * [O; I_k]
    //       = (Q1 ... Qn) * [O; I_k]
    //
    
    public var kernelMatrix: MatrixIF<Impl, m, anySize>  {
        let (r, m) = (rank, size.cols)
        if r == m {
            return .zero(size: (m, 0))
        }
        
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return .colUnits(
                size: (m, m - r),
                indices: (r ..< m)
            )
            .appliedRowOperations(
                colOps.reversed().map{ $0.asRowOperation }
            )
            
        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // Returns a matrix C consisting of the basis vectors of the free part of Coker(A).
    //
    //     P * A * Q = [ L | O   ]
    //                 [ * | O_k ]
    //
    //  =>  C_f = P^-1 * S,
    //      where S = (ei1 ... eil) consists of unit col-vectors that are complementary w.r.t. L.
    
    public var freeCokernelMatrix: MatrixIF<Impl, n, anySize> {
        let (r, n) = (rank, size.rows)
        if r == n {
            return .zero(size: (n, 0))
        }
        
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return .colUnits(
                size: (n, n - r),
                indices: complementRowsOfL
            )
            .appliedRowOperations(rowOpsInverse)
            
        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    fileprivate var complementRowsOfL: [Int] {
        let n = size.rows
        let occupied = Set(headEntries.map{ $0.row })
        return (0 ..< n).subtract(occupied)
    }
}

// MARK: Solutions

extension MatrixEliminationResult where R: EuclideanRing {
    // Find a solution x to: Ax = b.
    // With PAQ = B,
    //
    //    Ax = b  <==>  (PAQ) Q^{-1}x = Pb
    //            <==>    B      y    = Pb
    //
    // where y = Q^{-1}x <==> x = Qy.
    
    public func solve(_ b: MatrixIF<Impl, n, _1>) -> MatrixIF<Impl, m, _1>? {
        assert(isDiagonal) // TODO support non diagonal cases.
        
        let n = size.rows
        let r = rank
        let Pb = b.appliedRowOperations(rowOps)
        let diag = headEntries.map{ $0.value }
        
        if !Pb[r ..< n].isZero {
            return nil
        }
        
        if diag.enumerated().contains(where: { (i, d) in
            !d.divides(Pb[i])
        }) {
            return nil
        }
        
        let y = MatrixIF<Impl, m, _1>(size: (size.cols, 1)) { setEntry in
            diag.enumerated().forEach { (i, d) in
                setEntry(i, 0, Pb[i] / d)
            }
        }
        
        // return Qy
        return y.appliedRowOperations(
            colOps.reversed().map{ $0.asRowOperation }
        )
    }
    
    // Given z \in Span(Z), solve Zx = z.
    // We have
    //
    //   Zx = (Q1 ... Qn) [0; x],
    //
    // so it suffices to solve
    //
    //   [0; x] = (Qn^-1 ... Q1^-1) z.
    //
    public func solveKernel(_ z: ColVector<R, m>) -> AnySizeVector<R>? {
        assert(z.size.rows == size.cols)
        
        let (r, m) = (rank, size.cols)
        let w = z.appliedRowOperations(colOps.map{ $0.inverse.asRowOperation })
        if w[0 ..< r].isZero {
            return w[r ..< m]
        } else {
            return nil
        }
    }
    
    // Given z \in Span(C_f), solve C_f x = z.
    // We have
    //
    //   C_f x = P^-1 * Sx,
    //
    // so it suffices to solve
    //
    //   Sx = Σ_j x_j e_ij = Pz.
    //
    public func solveFreeCokernel(_ z: ColVector<R, n>) -> AnySizeVector<R>? {
        assert(z.size.rows == size.rows)
        
        let indices = complementRowsOfL
        let l = indices.count // n - rank
        
        let w = z.appliedRowOperations(rowOps)
        if w.nonZeroColEntries.map( {$0.row} ).subtract(Set(indices)).isEmpty {
            return .init(
                size: l,
                colEntries: indices.enumerated().map{ (j, i) in (j, w[i]) }
            )
        } else {
            return nil
        }
    }
}

// MARK: SNF specific

extension MatrixEliminationResult where R: EuclideanRing {
    public var divisors: [R] {
        guard form == .Smith else {
            fatalError("unavailable")
        }
        
        return headEntries.map{ $0.value }
    }

    public var nonUnitDivisors: [R] {
        guard form == .Smith else {
            fatalError("unavailable")
        }
        
        return divisors.exclude{ $0.isIdentity }
    }
    
    public var torCokernelMatrix: MatrixIF<Impl, n, anySize> {
        guard form == .Smith else {
            fatalError("unavailable")
        }
        
        let n = size.rows
        let r = rank
        let l = nonUnitDivisors.count

        if l == 0 {
            return .zero(size: (n, 0))
        }
        
        return .colUnits(
            size: (n, l),
            indices: (r - l ..< r)
        )
        .appliedRowOperations(rowOpsInverse)
    }
    
    // Given z \in Span(C_t), solve C_t x = z.
    // We have
    //
    //   C_t x = P^-1 * [0; x; 0],
    //
    // so it suffices to solve
    //
    //   [0; x; 0] = Pz.
    //
    public func solveTorCokernel(_ z: ColVector<R, n>) -> ColVector<R, anySize>? {
        guard form == .Smith else {
            fatalError("unavailable")
        }
        
        assert(z.size.rows == size.rows)
        
        let n = size.rows
        let r = rank
        let s = r - nonUnitDivisors.count
        
        let w = z.appliedRowOperations(rowOps)
        if w[0 ..< s].isZero && w[r ..< n].isZero {
            let d = divisors.exclude{ $0.isIdentity }
            return w[s ..< r].mapNonZeroEntries {
                (i, _, a) in a % d[i]
            }
        } else {
            return nil
        }
    }
}
