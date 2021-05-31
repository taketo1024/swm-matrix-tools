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
    public func permuteByPivots() -> (MatrixIF<Impl, n, m>, Permutation<n>, Permutation<m>) {
        let pf = MatrixPivotFinder(self)
        let pivots = pf.findPivots()
        
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
    public func eliminate(form: MatrixEliminationForm = .Diagonal, preprocess: Bool = false) -> MatrixEliminationResult<Impl, n, m> {
        if preprocess {
            let (B, p, q) = permuteByPivots()
            let e = B.eliminate(form: form, preprocess: false)
            return e.precompose(
                rowOps: p.asRowOps(),
                colOps: q.asColOps()
            )
        }
        
        let (type, transpose) = eliminatorType(form)
        let data = !transpose
            ? MatrixEliminationData(self)
            : MatrixEliminationData(self.transposed) // TODO directly pass tranposed entries
        
        let e = type.init(data: data)
        e.run()
        
        return !transpose
            ? e.result(as: MatrixEliminationResult.self)
            : e.result(as: MatrixEliminationResult.self).transposed
    }
    
    private func eliminatorType(_ form: MatrixEliminationForm) -> (MatrixEliminator<BaseRing>.Type, Bool) {
        switch form {
        case .RowEchelon:
            return (RowEchelonEliminator.self, false)
        case .ColEchelon:
            return (RowEchelonEliminator.self, true)
        case .RowHermite:
            return (ReducedRowEchelonEliminator.self, false)
        case .ColHermite:
            return (ReducedRowEchelonEliminator.self, true)
        case .Diagonal:
            return (DiagonalEliminator.self, false)
        case .Smith:
            return (SmithEliminator.self, false)
        default:
            return (MatrixEliminator.self, false)
        }
    }
}

extension MatrixEliminator {
    public convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>) where Impl.BaseRing == R {
        self.init(data: MatrixEliminationData(A))
    }
}

extension MatrixPivotFinder {
    public convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>) where Impl.BaseRing == R {
        self.init(data: MatrixEliminationData(A))
    }
}

extension MatrixEliminationData {
    convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>) where Impl.BaseRing == R {
        self.init(size: A.size, entries: A.nonZeroEntries)
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
