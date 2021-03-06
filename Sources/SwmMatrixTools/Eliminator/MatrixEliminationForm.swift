//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/25.
//

import SwmCore

public enum MatrixEliminationForm {
    case none
    case RowEchelon
    case ColEchelon
    case Diagonal
    case Smith
    
    var transposed: Self {
        switch self {
        case .Diagonal, .Smith, .none:
            return self
        case .RowEchelon:
            return .ColEchelon
        case .ColEchelon:
            return .RowEchelon
        }
    }
}
