//
//  LUFactorizable.swift
//  
//
//  Created by Taketo Sano on 2021/06/01.
//

import SwmCore

public protocol LUFactorizable: MatrixImpl {
    // solve L * x = b.
    static func solveLowerTriangular(_ L: Self, _ b: Self) -> Self?

    // solve U * x = b.
    static func solveUpperTriangular(_ U: Self, _ b: Self) -> Self?
}
