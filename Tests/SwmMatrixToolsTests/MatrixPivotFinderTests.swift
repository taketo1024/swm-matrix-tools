//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/28.
//

import XCTest
import SwmCore
@testable import SwmMatrixTools

class MatrixPivotFinderTests: XCTestCase {
    
    typealias M<n: SizeType, m: SizeType> = Matrix<Int, n, m>
    
    func testPivotFinder() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let pf = MatrixPivotFinder(A)
        XCTAssertEqual(pf.pivots.count, 5)
    }
    
    func testPermuteByPivots() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let pivs = A.findPivots()
        let (B, p, q) = A.permute(byPivots: pivs)
        
        XCTAssertTrue(pivs.count >= 4)
        XCTAssertEqual(A.permute(rowsBy: p, colsBy: q), B)
        XCTAssertTrue(
            B.submatrix(rowRange: 0 ..< pivs.count, colRange: 0 ..< pivs.count)
                .nonZeroEntries
                .allSatisfy{ (i, j, a) in i < j || (i == j && a.isIdentity) }
        )
    }
}
