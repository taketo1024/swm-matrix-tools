//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/26.
//

import SwmCore

extension MatrixIF where BaseRing: EuclideanRing {
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
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
