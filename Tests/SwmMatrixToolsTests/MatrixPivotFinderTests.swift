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
        let pivs = MatrixPivotFinder(A).findPivots()
        XCTAssertEqual(pivs.count, 5)
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
        let (B, p, q) = A.permuteByPivots()
        
        XCTAssertEqual(A.permute(rowsBy: p, colsBy: q), B)
        
        // check upper-triangular
        XCTAssertTrue(
            B.submatrix(rowRange: 0 ..< 5, colRange: 0 ..< 5)
                .nonZeroEntries
                .allSatisfy{ (i, j, _) in i <= j }
        )
    }
}
