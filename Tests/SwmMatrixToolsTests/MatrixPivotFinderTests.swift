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
    
    func testSample() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let pivs = A.findPivots()
        
        XCTAssertEqual(pivs.count, 5)
        
        let (p, q) = A.permutations(forPivots: pivs)
        let B = A.permute(rowsBy: p, colsBy: q)
        
        XCTAssertTrue(B.submatrix(rowRange: 0 ..< 5, colRange: 0 ..< 5).nonZeroEntries.allSatisfy{ (i, j, _) in i <= j }) // upper-triangular
    }
}
