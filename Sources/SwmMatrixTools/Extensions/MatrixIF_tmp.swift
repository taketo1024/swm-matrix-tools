//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/26.
//

import SwmCore

// TODO remove these. too much complicated.
extension MatrixEliminationResult {
    
    public func left(restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, anySize, n> {
        composeRowOps(rowOps, restrictedToRows: rowRange)
    }
    
    public func leftInverse(restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, anySize> {
        composeRowOps(rowOpsInverse, restrictedToCols: colRange)
    }
    
    public func right(restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, m, anySize> {
        composeColOps(colOps, restrictedToCols: colRange)
    }
    
    public func rightInverse(restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, anySize, m> {
        composeColOps(colOpsInverse, restrictedToRows: rowRange)
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
    
    //  Given row ops [P1, ..., Pn],
    //  produce P = (Pn ... P1) * I by applying the row-ops from P1 to Pn.
    
    private func composeRowOps<n, m, S>(_ ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        composeRowOps(size: size.rows, ops: ops, restrictedToCols: colRange)
    }
    
    //  Given row ops [P1, ..., Pn],
    //  produce P = I * (Pn ... P1) by applying the corresponding col-ops from Pn to P1.
    //  Compute by P^t = (P1^t ... Pn^t) * I^t,
    
    private func composeRowOps<n, m, S>(_ ops: S, restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        composeRowOps(size: size.rows, ops: ops.reversed().map{ $0.transposed.asRowOperation }, restrictedToCols: rowRange).transposed
    }
    
    //  Given col ops [Q1, ..., Qn],
    //  produce Q = I * (Q1 ... Qn) by applying the col-ops from Q1 to Qn.
    //  Compute by Q^t = (Qn^t ... Q1^t) * I^t.

    private func composeColOps<n, m, S>(_ ops: S, restrictedToRows rowRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == ColElementaryOperation<R> {
        composeRowOps(size: size.cols, ops: ops.map{ $0.transposed }, restrictedToCols: rowRange).transposed
    }
    
    //  Given col ops [Q1, ..., Qn],
    //  produce Q = (Q1 ... Qn) * I by applying the corresponding row-ops from Qn to Q1.
    
    private func composeColOps<n, m, S>(_ ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == ColElementaryOperation<R> {
        composeRowOps(size: size.cols, ops: ops.reversed().map{ $0.asRowOperation }, restrictedToCols: colRange)
    }
    
    private func composeRowOps<n, m, S>(size n: Int, ops: S, restrictedToCols colRange: Range<Int>) -> MatrixIF<Impl, n, m> where S: Sequence, S.Element == RowElementaryOperation<R> {
        .init(
            size: (n, colRange.endIndex - colRange.startIndex),
            entries: colRange.map { i in (i, i - colRange.startIndex, .identity) }
        )
        .appliedRowOperations(ops)
    }
}
