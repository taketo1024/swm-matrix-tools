//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/10.
//

import SwmCore

public protocol ComputationalRing: Ring {
    var computationalWeight: Double { get } // used for matrix elimination
}

extension ComputationalRing where Self: FractionField, Base: ComputationalRing {
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : max(numerator.computationalWeight, denominator.computationalWeight)
    }
}

extension Int: ComputationalRing {
    @inlinable
    public var computationalWeight: Double {
        Double(Swift.abs(self))
    }
}

extension RationalNumber: ComputationalRing {}

extension RealNumber: ComputationalRing {
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double( Swift.abs( max(self, 1/self) ) )
    }
}

extension 𝐅₂: ComputationalRing {
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : 1
    }
}

extension Polynomial: ComputationalRing where BaseRing: ComputationalRing {
    @inlinable
    public var computationalWeight: Double {
        isZero ? 0 : Double(leadExponent + 1) * leadCoeff.computationalWeight
    }
}
