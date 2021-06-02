//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/26.
//

import SwmCore

extension MatrixIF {
    public var isLowerTriangular: Bool {
        impl.isLowerTriangular
    }
    
    public var isUpperTriangular: Bool {
        impl.isUpperTriangular
    }
    
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
    
    internal func appliedRowOperations<S>(_ ops: S) -> Self
    where S: Sequence, S.Element == RowElementaryOperation<BaseRing> {
        MatrixEliminationData(self).applyAll(ops).resultAs(Self.self)
    }
    
    internal func appliedColOperations<S>(_ ops: S) -> Self
    where S: Sequence, S.Element == ColElementaryOperation<BaseRing> {
        transposed.appliedRowOperations(ops.map{ $0.transposed }).transposed
    }
}

extension MatrixIF {
    public func findPivots(mode: PivotMode = .rowBased) -> (pivots: [(Int, Int)], P: Permutation<n>, Q: Permutation<m>) {
        let pf = MatrixPivotFinder(self, mode: mode)
        pf.run()
        
        return (
            pivots: pf.pivots,
            P: pf.rowPermutation.as(Permutation.self),
            Q: pf.colPermutation.as(Permutation.self)
        )
    }
}

extension MatrixIF where BaseRing: EuclideanRing {
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
        MatrixEliminator.eliminate(self, form: form)
    }
}

extension MatrixIF where Impl: LUFactorizable {
    public func LUfactorize() -> LUFactorizationResult<Impl, n, m> {
        let (P, Q, L, U) = impl.LUfactorize()
        return LUFactorizationResult(
            P: P.as(Permutation.self),
            Q: Q.as(Permutation.self),
            L: .init(L),
            U: .init(U)
        )
    }
    
    public static func solveLowerTrapezoidal<k>(_ L: Self, _ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k>? {
        Impl.solveLowerTrapezoidal(L.impl, b.impl).flatMap{ .init($0) }
    }
    
    public static func solveUpperTrapezoidal<k>(_ L: Self, _ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k> {
        .init(Impl.solveUpperTrapezoidal(L.impl, b.impl))
    }
}

extension MatrixIF where Impl: LUFactorizable, n == m {
    public static func solveLowerTriangular<k>(_ L: Self, _ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k> {
        .init(Impl.solveLowerTriangular(L.impl, b.impl))
    }
    
    public static func solveUpperTriangular<k>(_ U: Self, _ b: MatrixIF<Impl, n, k>) -> MatrixIF<Impl, m, k> {
        .init(Impl.solveUpperTriangular(U.impl, b.impl))
    }
}
