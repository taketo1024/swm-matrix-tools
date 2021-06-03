//
//  LinkedList.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/27.
//

internal final class LinkedList<Element>: Sequence {
    fileprivate typealias RawPointer = UnsafeMutablePointer<Node>
    private var headPtr: RawPointer? = nil
    
    init<S: Sequence>(_ seq: S) where S.Element == Element {
        var head: RawPointer?
        var prev: RawPointer?
        
        for e in seq {
            let p = RawPointer.new( Node(element: e) )

            if head == nil {
                head = p
            }
            
            prev?.pointee.next = p
            prev = p
        }
        
        self.headPtr = head
    }
    
    convenience init() {
        self.init([])
    }
    
    deinit {
        removeAll()
    }
    
    var isEmpty: Bool {
        headPtr == nil
    }
    
    var head: NodePointer? {
        headPtr.flatMap{ .init($0) }
    }
    
    func insertHead(_ element: Element) {
        headPtr = RawPointer.new( Node(element: element, next: headPtr) )
    }
    
    func dropHead() {
        guard let head = headPtr else {
            return
        }
        defer { head.delete() }
        self.headPtr = head.pointee.next
    }
    
    func removeAll() {
        var p = headPtr
        while p != nil {
            let next = p?.pointee.next
            
            p!.delete()
            p = next
        }
        headPtr = nil
    }
    
    func modifyEach(_ map: (inout Element) -> Void) {
        var current = head
        while var p = current {
            map(&(p.element))
            current = p.next
        }
    }
    
    func first(where hit: (Element) -> Bool, while cond: (Element) -> Bool) -> Element? {
        var current = head
        while true {
            guard let p = current, cond(p.element) else {
                break
            }
            if hit(p.element) {
                return p.element
            }
            current = p.next
        }
        return nil
    }
    
    fileprivate struct Node {
        var element: Element
        var next: RawPointer? = nil
        
        mutating func insertNext(_ e: Element) {
            let p = RawPointer.new( Node(element: e, next: self.next) )
            self.next = p
        }
        
        mutating func dropNext() {
            guard let drop = next else {
                return
            }
            self.next = drop.pointee.next
            drop.delete()
        }
    }
    
    func makeIterator() -> ElementIterator {
        ElementIterator(headPtr?.pointee)
    }
    
    struct NodePointer: Equatable {
        private var ptr: RawPointer
        fileprivate init(_ ptr: RawPointer) {
            self.ptr = ptr
        }
        
        var element: Element {
            get {
                ptr.pointee.element
            } set {
                ptr.pointee.element = newValue
            }
        }
        
        var hasNext: Bool {
            ptr.pointee.next != nil
        }
        
        var next: NodePointer? {
            ptr.pointee.next.flatMap{ .init($0) }
        }
        
        func insertNext(_ e: Element) {
            ptr.pointee.insertNext(e)
        }
        
        func dropNext() {
            ptr.pointee.dropNext()
        }
    }
    
    struct ElementIterator: IteratorProtocol {
        private var current: Node?
        fileprivate init(_ start: Node?) {
            current = start
        }
        
        mutating func next() -> Element? {
            defer { current = current?.next?.pointee }
            return current?.element
        }
    }
}

private extension UnsafeMutablePointer {
    static func new(_ entity: Pointee) -> Self {
        let p = allocate(capacity: 1)
        p.initialize(to: entity)
        return p
    }
    
    func delete() {
        self.deinitialize(count: 1)
        self.deallocate()
    }
}
