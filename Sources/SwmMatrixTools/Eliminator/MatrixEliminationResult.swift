//
//  MatrixEliminationResult.swift
//  Sample
//
//  Created by Taketo Sano on 2018/04/26.
//

import SwmCore

public struct MatrixEliminationResult<Impl, n, m>
where Impl: MatrixImpl, Impl.BaseRing: ComputationalRing, n: SizeType, m: SizeType {
    public typealias R = Impl.BaseRing
    
    public let form: MatrixEliminationForm
    public let result: MatrixIF<Impl, n, m>
    public let headEntries: [MatrixEntry<R>]
    internal let rowOps: [RowElementaryOperation<R>]
    internal let colOps: [ColElementaryOperation<R>]
    
    private let cache: Cache<String, MatrixIF<Impl, anySize, anySize>> = .empty
    
    internal init(form: MatrixEliminationForm, result: MatrixIF<Impl, n, m>, headEntries: [MatrixEntry<R>], rowOps: [RowElementaryOperation<R>], colOps: [ColElementaryOperation<R>]) {
        self.form = form
        self.result = result
        self.headEntries = headEntries
        self.rowOps = rowOps
        self.colOps = colOps
    }
    
    public var size: MatrixSize {
        result.size
    }
    
    public var isSquare: Bool {
        size.rows == size.cols
    }
    
    public var rank: Int {
        assert(form != .none)
        return headEntries.count
    }
    
    public var nullity: Int {
        size.cols - rank
    }
    
    public var left: MatrixIF<Impl, n, n> {
        cache.getOrSet(key: "left") {
            AnySizeMatrixIF
                .identity(size: (size.rows, size.rows))
                .appliedRowOperations(rowOps)
        }.as(MatrixIF.self)
    }
    
    public var leftInverse: MatrixIF<Impl, n, n> {
        cache.getOrSet(key: "leftInverse") {
            AnySizeMatrixIF
                .identity(size: (size.rows, size.rows))
                .appliedRowOperations(rowOpsInverse)
        }.as(MatrixIF.self)
    }
    
    public var right: MatrixIF<Impl, m, m> {
        cache.getOrSet(key: "right") {
            AnySizeMatrixIF
                .identity(size: (size.cols, size.cols))
                .appliedColOperations(colOps)
        }.as(MatrixIF.self)
    }
    
    public var rightInverse: MatrixIF<Impl, m, m> {
        cache.getOrSet(key: "rightInverse") {
            AnySizeMatrixIF
                .identity(size: (size.cols, size.cols))
                .appliedColOperations(colOpsInverse)
        }.as(MatrixIF.self)
    }
    
    private var rowOpsInverse: [RowElementaryOperation<R>] {
        rowOps.reversed().map{ $0.inverse }
    }

    private var colOpsInverse: [ColElementaryOperation<R>] {
        colOps.reversed().map{ $0.inverse }
    }
    
    public var transposed: MatrixEliminationResult<Impl, m, n> {
        return .init(
            form: form.transposed,
            result: result.transposed,
            headEntries: headEntries.map{ (i, j, a) in (j, i, a) },
            rowOps: colOps.map{ $0.transposed },
            colOps: rowOps.map{ $0.transposed }
        )
    }
    
    public func precompose(rowOps: [RowElementaryOperation<R>], colOps: [ColElementaryOperation<R>]) -> Self {
        .init(
            form: form,
            result: result,
            headEntries: headEntries,
            rowOps: rowOps + self.rowOps,
            colOps: colOps + self.colOps
        )
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
        return (result.isIdentity) ? right * left : nil
    }
}

extension MatrixEliminationResult where R: EuclideanRing {
    // eliminate again
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
        result.eliminate(form: form)
            .precompose(rowOps: rowOps, colOps: colOps)
    }
}

// MARK: Associated matrices

extension MatrixEliminationResult {
    
    // Returns a matrix consisting of a basis (col-vectors) of Im(A).
    
    public var imageMatrix: MatrixIF<Impl, n, anySize> {
        switch form {
        case .Diagonal, .Smith, .ColEchelon:
            return imageMatrix_colEchelon
            
        case .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // With PAQ = [L|O], Im(PA) is spanned by L.
    // => B = P^-1 * L.
    
    private var imageMatrix_colEchelon: MatrixIF<Impl, n, anySize> {
        let A = result // must obtain outside cache-sync
        let Pinv = leftInverse
        
        return cache.getOrSet(key: "image") {
            let (r, n) = (rank, size.rows)
            if r == 0 {
                return .zero(size: (n, 0))
            }
            
            let L = A.submatrix(colRange: 0 ..< r)
            return (Pinv * L).asAnySizeMatrix
            
        }.as(MatrixIF.self)
    }
    
    // Returns a matrix Z consisting of a basis vectors of Ker(A).
    
    public var kernelMatrix: MatrixIF<Impl, m, anySize>  {
        switch form {
        case .Diagonal, .Smith, .ColEchelon:
            return kernelMatrix_colEchelon
            
        case .RowEchelon:
            fatalError("not supported yet.")
            
        default:
            fatalError("unavailable.")
        }
    }
    
    // With PAQ = [L|O], Ker(AQ) is spanned by [O|I].
    // => Ker(A) = Q * [O|I].
    
    private var kernelMatrix_colEchelon: MatrixIF<Impl, m, anySize>  {
        let Q = right // must obtain outside cache-sync
        
        return cache.getOrSet(key: "kernel") {
            let (r, m) = (rank, size.cols)
            if r == m {
                return .zero(size: (m, 0))
            }
            
            let I = MatrixIF<Impl, m, anySize>.colUnits(
                size: (m, m - r),
                indices: (r ..< m)
            )
            return (Q * I).asAnySizeMatrix
            
        }.as(MatrixIF.self)
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
    
    public func solve(_ b: ColVectorIF<Impl, n>) -> ColVectorIF<Impl, m>? {
        assert(form == .Diagonal) // TODO support non diagonal cases.
        
        let n = size.rows
        let r = rank
        
        let P = left
        let Pb = P * b
        let diag = headEntries.map{ $0.value }
        
        if !Pb[r ..< n].isZero {
            return nil
        }
        
        if diag.enumerated().contains(where: { (i, d) in
            !d.divides(Pb[i])
        }) {
            return nil
        }
        
        let Q = right
        let y = ColVectorIF<Impl, m>(size: (size.cols, 1)) { setEntry in
            diag.enumerated().forEach { (i, d) in
                setEntry(i, 0, Pb[i] / d)
            }
        }
        
        return Q * y
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
    public func solveKernel(_ z: ColVectorIF<Impl, m>) -> ColVectorIF<Impl, anySize>? {
        assert(z.size.rows == size.cols)
        
        let (r, m) = (rank, size.cols)
        let Qinv = rightInverse
        let w = Qinv * z
        
        if w[0 ..< r].isZero {
            return w[r ..< m]
        } else {
            return nil
        }
    }
}
