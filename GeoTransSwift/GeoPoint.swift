//
//  GeoPoint.swift
//  GeoTransSwift
//
//  Created by steve.ham on 2022/01/20.
//

import Foundation

class GeoPoint {
    var longitude: Double
    var latitude: Double
    var altitude: Double
    
    init() {
        longitude = 0.0
        latitude = 0.0
        altitude = 0.0
    }
    
    init(longitude: Double, latitude: Double) {
        self.longitude = longitude
        self.latitude = latitude
        altitude = 0.0
    }
    
    init(longitude: Double, latitude: Double, altitude: Double) {
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude
    }
}
