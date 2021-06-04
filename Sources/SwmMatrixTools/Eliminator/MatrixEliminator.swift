//
//  MatrixEliminator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/06/09.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import SwmCore

public class MatrixEliminator<R: Ring> {
    internal let data: MatrixEliminationData<R>
    internal var rowOps: [RowElementaryOperation<R>]
    internal var colOps: [ColElementaryOperation<R>]
    
    public private(set) var aborted: Bool = false
    public var debug: Bool = false

    required init(data: MatrixEliminationData<R>) {
        self.data = data
        self.rowOps = []
        self.colOps = []
    }
    
    public final func run() {
        log("Start: \(self)")
        
        prepare()
        
        var itr = 0
        while !aborted && !isDone() {
            log("\(self) iteration: \(itr)")
            
            if debug {
                printCurrentMatrix()
            }
            
            iteration()
            
            log("")
            itr += 1
        }
        
        finalize()
        
        log("Done:  \(self), \(rowOps.count + colOps.count) steps")
        
        if debug {
            printCurrentMatrix()
        }
    }
    
    public func result<Impl, n, m>(as: MatrixEliminationResult<Impl, n, m>.Type) -> MatrixEliminationResult<Impl, n, m> where Impl.BaseRing == R {
        .init(
            form: form,
            result: data.resultAs(MatrixIF.self),
            headEntries: data.headEntries,
            rowOps: rowOps,
            colOps: colOps
        )
    }
    
    public var description: String {
        "\(type(of: self))"
    }
    
    // MARK: Internal methods

    final func subrun(_ type: MatrixEliminator<R>.Type, transpose: Bool = false) {
        let e = type.init(data: data)
        if transpose {
            e.transpose()
        }
        
        e.run()
        
        if transpose {
            e.transpose()
        }
        
        rowOps += e.rowOps
        colOps += e.colOps
    }
    
    final func abort() {
        log("Aborted.")
        aborted = true
    }
    
    final func apply(_ s: RowElementaryOperation<R>) {
        data.apply(s)
        append(s)
    }

    final func append(_ s: RowElementaryOperation<R>) {
        rowOps.append(s)
        log("\(s)")
    }
    
    final func append(_ s: ColElementaryOperation<R>) {
        colOps.append(s)
        log("\(s)")
    }
    
    final func append(_ s: [RowElementaryOperation<R>]) {
        rowOps.append(contentsOf: s)
        log(s.map{ "\($0)"}.joined(separator: "\n"))
    }
    
    final func transpose() {
        log("Transpose: \(self)")
        data.transpose()
        (rowOps, colOps) = (colOps.map{ $0.transposed }, rowOps.map{ $0.transposed })
    }
    
    final func log(_ msg: @autoclosure () -> String) {
        if debug {
            print(msg())
        }
    }
    
    final func printCurrentMatrix() {
        if data.size.rows > 100 || data.size.cols > 100 {
            return
        }
        print("\n", data.resultAs(AnySizeMatrix.self).detailDescription, "\n")
    }

    // MARK: Methods to be overridden
    
    public var form: MatrixEliminationForm { .none }
    func prepare() {}
    func isDone() -> Bool { true }
    func iteration() {}
    func finalize() {}
}

extension MatrixEliminator where R: EuclideanRing {
    public static func eliminate<Impl, n, m>(_ A: MatrixIF<Impl, n, m>, form: MatrixEliminationForm) -> MatrixEliminationResult<Impl, n, m>
    where Impl: MatrixImpl, Impl.BaseRing == R { 
        let (type, transpose) = eliminatorType(form)
        let data = MatrixEliminationData(A, transpose: transpose)
        let e = type.init(data: data)
        
        e.run()
        
        return !transpose
            ? e.result(as: MatrixEliminationResult.self)
            : e.result(as: MatrixEliminationResult.self).transposed
    }
    
    private static func eliminatorType(_ form: MatrixEliminationForm) -> (type: MatrixEliminator<R>.Type, transpose: Bool) {
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
