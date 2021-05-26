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
    
    public var diagonalEntries: [R] {
        headEntries.map{ $0.value }
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
                * diagonalEntries.multiply()
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

extension MatrixEliminationResult {
    
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
            return .rowUnits(
                size: (m, m - r),
                indices: (0 ..< m - r)
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
}

extension MatrixEliminationResult where R: EuclideanRing {
    // Find a solution x to: Ax = b.
    // With PAQ = B,
    //
    //    Ax = b  <==>  (PAQ) Q^{-1}x = Pb
    //            <==>    B      y    = Pb
    //
    // where y = Q^{-1}x <==> x = Qy.
    
    public func invert(_ b: MatrixIF<Impl, n, _1>) -> MatrixIF<Impl, m, _1>? {
        assert(isDiagonal)
        let r = rank
        let P = left
        let Pb = P * b
        
        if diagonalEntries.enumerated().contains(where: { (i, d) in
            (d.isZero && !Pb[i].isZero) || (!d.isZero && !Pb[i].isDivible(by: d))
        }) {
            return nil // no solution
        }
        
        if Pb.nonZeroEntries.contains(where: { (i, _, a) in i >= r && !a.isZero } ) {
            return nil // no solution
        }
        
        let Q = right
        let y = MatrixIF<Impl, m, _1>(size: (size.cols, 1)) { setEntry in
            diagonalEntries.enumerated().forEach { (i, d) in
                if !d.isZero {
                    setEntry(i, 0, Pb[i] / d)
                }
            }
        }
        return Q * y
    }
}
