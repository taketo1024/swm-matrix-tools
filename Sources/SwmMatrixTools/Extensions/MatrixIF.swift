//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/26.
//

import SwmCore

extension MatrixIF {
    public static func rowUnits<S>(size: MatrixSize, indices: S) -> Self
    where S: Sequence, S.Element == Int {
        assert(size.rows == indices.count)
        return .init(
            size: size,
            entries: indices.enumerated().map{ (i, idx) in (i, idx, .identity) }
        )
    }
    
    public static func colUnits<S>(size: MatrixSize, indices: S) -> Self
    where S: Sequence, S.Element == Int {
        assert(size.cols == indices.count)
        return .init(
            size: size,
            entries: indices.enumerated().map{ (j, idx) in (idx, j, .identity) }
        )
    }
    
    public func appliedRowOperations<S>(_ ops: S) -> Self
    where S: Sequence, S.Element == RowElementaryOperation<BaseRing> {
        MatrixEliminationData(
            size: size,
            entries: nonZeroEntries
        ).applyAll(ops).resultAs(Self.self)
    }
    
    public func appliedColOperations<S>(_ ops: S) -> Self
    where S: Sequence, S.Element == ColElementaryOperation<BaseRing> {
        transposed.appliedRowOperations(ops.map{ $0.transposed }).transposed
    }
}

extension MatrixIF {
    public func findPivots() -> [(Int, Int)] {
        MatrixPivotFinder(self).pivots
    }
    
    public func permute(byPivots pivots: [(Int, Int)]) -> (MatrixIF<Impl, n, m>, Permutation<n>, Permutation<m>) {
        let p: Permutation<n> = asPermutation(size.rows, pivots.map{ $0.0 })
        let q: Permutation<m> = asPermutation(size.cols, pivots.map{ $0.1 })
        
        return (permute(rowsBy: p, colsBy: q), p, q)
    }
    
    private func asPermutation<n>(_ length: Int, _ order: [Int]) -> Permutation<n> {
        let remain = Set(0 ..< length).subtracting(order)
        let p = Permutation<n>(length: length, indices: order + remain.sorted())
        return p.inverse!
    }
}

extension MatrixIF where BaseRing: EuclideanRing {
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
        MatrixEliminator.eliminate(self, form: form)
    }
}

extension Permutation {
    func asRowOps<R: Ring>() -> [RowElementaryOperation<R>] {
        transpositionDecomposition.map { (i, j) in .SwapRows(i, j) }
    }

    func asColOps<R: Ring>() -> [ColElementaryOperation<R>] {
        transpositionDecomposition.map { (i, j) in .SwapCols(i, j) }
    }
}
