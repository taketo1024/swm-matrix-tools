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

    let rowOps: [RowElementaryOperation<R>]
    let colOps: [ColElementaryOperation<R>]
    
    internal init(form: MatrixEliminationForm, size: MatrixSize, entries: AnySequence<MatrixEntry<R>>, headEntries: AnySequence<MatrixEntry<R>>, rowOps: [RowElementaryOperation<R>], colOps: [ColElementaryOperation<R>]) {
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
    
    // returns P of: P * A * Q = B

    public var left: MatrixIF<Impl, n, n> {
        composeRowOps(rowOps, restrictedToRows: 0 ..< size.rows)
    }
    
    public func left(restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, anySize, n> {
        composeRowOps(rowOps, restrictedToRows: rowRange)
    }
    
    // returns P^{-1} of: P * A * Q = B
    
    public var leftInverse: MatrixIF<Impl, n, n> {
        composeRowOps(rowOpsInverse, restrictedToCols: 0 ..< size.rows)
    }
    
    public func leftInverse(restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, anySize> {
        composeRowOps(rowOpsInverse, restrictedToCols: colRange)
    }
    
    // returns Q of: P * A * Q = B
    
    public var right: MatrixIF<Impl, m, m> {
        composeColOps(colOps, restrictedToRows: 0 ..< size.cols)
    }
    
    public func right(restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, m, anySize> {
        composeColOps(colOps, restrictedToCols: colRange)
    }
    
    // returns Q^{-1} of: P * A * Q = B
    
    public var rightInverse: MatrixIF<Impl, m, m> {
        composeColOps(colOpsInverse, restrictedToRows: 0 ..< size.cols)
    }
    
    public func rightInverse(restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, anySize, m> {
        composeColOps(colOpsInverse, restrictedToRows: rowRange)
    }
    
    // Returns the matrix Z consisting of the basis vectors of Ker(A).
    // With k = m - rank,
    //
    //     P * A * Q = [ L  | O   ]
    //                 [ L' | O_k ]
    //
    //  => Z = Q * [O; I_k]
    //       = Q[-, r ..< m]
    //
    
    public var kernelMatrix: MatrixIF<Impl, m, anySize>  {
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return right(restrictedToCols: rank ..< size.cols)
            
        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // Returns the transition matrix T of Z,
    // i.e. T * Z = I_k.
    //
    //     T = [O, I_k] Q^-1
    //       = Q^-1 [r ..< n; -]
    //
    // satisfies T * Z = [O, I_k] * [O; I_k] = I_k.
    
    public var kernelTransitionMatrix: MatrixIF<Impl, anySize, m> {
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return rightInverse(restrictedToRows: rank ..< size.cols)
            
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
    // In case L' = 0,
    //
    //     B = P^-1[-; 0 ..< r] * L.
    
    public var imageMatrix: MatrixIF<Impl, n, anySize> {
        let r = rank
        
        switch form {
        case .Diagonal, .Smith:
            return leftInverse(restrictedToCols: 0 ..< r) * result.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)

        case .ColHermite, .ColEchelon:
            return leftInverse * result.submatrix(colRange: 0 ..< r)

        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // Returns a matrix C consisting of the basis vectors of Coker(A).
    // With
    //
    //     P * A * Q = [ L | O   ]
    //                 [ * | O_k ]
    //
    // Coker(A) is spanned by [O; I_{n-r}].
    //
    //  => C = P^-1 * [O; I_{n-r}]
    //       = P^-1 [r ..< n, -]

    public var cokernelMatrix: MatrixIF<Impl, n, anySize> {
        let (r, n) = (rank, size.rows)
        
        switch form {
        case .Diagonal, .Smith, .ColHermite, .ColEchelon:
            return leftInverse(restrictedToCols: r ..< n)

        case .RowHermite, .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    

    var rowOpsInverse: [RowElementaryOperation<R>] {
        rowOps.reversed().map{ $0.inverse }
    }

    var colOpsInverse: [ColElementaryOperation<R>] {
        colOps.reversed().map{ $0.inverse }
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

private extension MatrixEliminationResult {
    
    //  Given row ops [P1, ..., Pn],
    //  produce P = (Pn ... P1) * I by applying the row-ops from P1 to Pn.
    
    func composeRowOps<n, m, S>(_ ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        composeRowOps(size: size.rows, ops: ops, restrictedToCols: colRange)
    }
    
    //  Given row ops [P1, ..., Pn],
    //  produce P = I * (Pn ... P1) by applying the corresponding col-ops from Pn to P1.
    //  Compute by P^t = (P1^t ... Pn^t) * I^t,
    
    func composeRowOps<n, m, S>(_ ops: S, restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        composeRowOps(size: size.rows, ops: ops.reversed().map{ $0.transposed.asRowOperation }, restrictedToCols: rowRange).transposed
    }
    
    //  Given col ops [Q1, ..., Qn],
    //  produce Q = I * (Q1 ... Qn) by applying the col-ops from Q1 to Qn.
    //  Compute by Q^t = (Qn^t ... Q1^t) * I^t.

    func composeColOps<n, m, S>(_ ops: S, restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == ColElementaryOperation<R> {
        composeRowOps(size: size.cols, ops: ops.map{ $0.transposed }, restrictedToCols: rowRange).transposed
    }
    
    //  Given col ops [Q1, ..., Qn],
    //  produce Q = (Q1 ... Qn) * I by applying the corresponding row-ops from Qn to Q1.
    
    func composeColOps<n, m, S>(_ ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == ColElementaryOperation<R> {
        composeRowOps(size: size.cols, ops: ops.reversed().map{ $0.asRowOperation }, restrictedToCols: colRange)
    }
    
    private func composeRowOps<n, m, S>(size n: Int, ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        MatrixEliminationWorker(
            size: (n, colRange.endIndex - colRange.startIndex),
            entries: colRange.map { i in (i, i - colRange.startIndex, .identity) }
        ).applyAll(ops).resultAs(MatrixIF.self)
    }
}
