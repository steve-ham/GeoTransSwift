//
//  GeoTrans.swift
//  GeoTransSwift
//
//  Created by steve.ham on 2022/01/20.
//

import Foundation

class GeoTrans {
    
    static let shared = GeoTrans()
    
    enum CoordinateSystem: Int {
        case wgs84
        case katec
        case tm
        case grs80
        case utmk
    }
    
    private var m_Ind = Array(repeating: 0.0, count: 5)
    private var m_Es = Array(repeating: 0.0, count: 5)
    private var m_Esp = Array(repeating: 0.0, count: 5)
    private var src_m = Array(repeating: 0.0, count: 5)
    private var dst_m = Array(repeating: 0.0, count: 5)
    
    private let EPSLN = 0.0000000001
    private var m_arMajor = Array(repeating: 0.0, count: 5)
    private var m_arMinor = Array(repeating: 0.0, count: 5)
    
    private var m_arScaleFactor = Array(repeating: 0.0, count: 5)
    private var m_arLonCenter = Array(repeating: 0.0, count: 5)
    private var m_arLatCenter = Array(repeating: 0.0, count: 5)
    private var m_arFalseNorthing = Array(repeating: 0.0, count: 5)
    private var m_arFalseEasting = Array(repeating: 0.0, count: 5)
    
    private var datum_params = Array(repeating: 0.0, count: 3)
    
    private init() {
        m_arScaleFactor[CoordinateSystem.wgs84.rawValue] = 1
        m_arLonCenter[CoordinateSystem.wgs84.rawValue] = 0.0
        m_arLatCenter[CoordinateSystem.wgs84.rawValue] = 0.0
        m_arFalseNorthing[CoordinateSystem.wgs84.rawValue] = 0.0
        m_arFalseEasting[CoordinateSystem.wgs84.rawValue] = 0.0
        m_arMajor[CoordinateSystem.wgs84.rawValue] = 6378137.0
        m_arMinor[CoordinateSystem.wgs84.rawValue] = 6356752.3142
        
        m_arScaleFactor[CoordinateSystem.katec.rawValue] = 0.9999
        m_arLonCenter[CoordinateSystem.katec.rawValue] = 2.23402144255274 // 128
        m_arLatCenter[CoordinateSystem.katec.rawValue] = 0.663225115757845
        m_arFalseNorthing[CoordinateSystem.katec.rawValue] = 600000.0
        m_arFalseEasting[CoordinateSystem.katec.rawValue] = 400000.0
        m_arMajor[CoordinateSystem.katec.rawValue] = 6377397.155
        m_arMinor[CoordinateSystem.katec.rawValue] = 6356078.9633422494
        
        m_arScaleFactor[CoordinateSystem.tm.rawValue] = 1.0
        //this.m_arLonCenter[CoordinateSystem.tm.rawValue] = 2.21656815003280 // 127
        m_arLonCenter[CoordinateSystem.tm.rawValue] = 2.21661859489671 // 127.+10.485 minute
        m_arLatCenter[CoordinateSystem.tm.rawValue] = 0.663225115757845
        m_arFalseNorthing[CoordinateSystem.tm.rawValue] = 500000.0
        m_arFalseEasting[CoordinateSystem.tm.rawValue] = 200000.0
        m_arMajor[CoordinateSystem.tm.rawValue] = 6377397.155
        m_arMinor[CoordinateSystem.tm.rawValue] = 6356078.9633422494
        
        m_arScaleFactor[CoordinateSystem.grs80.rawValue] = 1.0//0.9999
        m_arLonCenter[CoordinateSystem.grs80.rawValue] = 2.21656815003280 // 127
        //m_arLonCenter[CoordinateSystem.grs80.rawValue] = 2.21661859489671 // 127.+10.485 minute
        m_arLatCenter[CoordinateSystem.grs80.rawValue] = 0.663225115757845
        m_arFalseNorthing[CoordinateSystem.grs80.rawValue] = 500000.0
        m_arFalseEasting[CoordinateSystem.grs80.rawValue] = 200000.0
        m_arMajor[CoordinateSystem.grs80.rawValue] = 6378137.0
        m_arMinor[CoordinateSystem.grs80.rawValue] = 6356752.3142
        
        m_arScaleFactor[CoordinateSystem.utmk.rawValue] = 0.9996//0.9999
        //m_arLonCenter[CoordinateSystem.utmk.rawValue] = 2.22534523630815 // 127.502890
        m_arLonCenter[CoordinateSystem.utmk.rawValue] = 2.22529479629277 // 127.5
        m_arLatCenter[CoordinateSystem.utmk.rawValue] = 0.663225115757845
        m_arFalseNorthing[CoordinateSystem.utmk.rawValue] = 2000000.0
        m_arFalseEasting[CoordinateSystem.utmk.rawValue] = 1000000.0
        m_arMajor[CoordinateSystem.utmk.rawValue] = 6378137.0
        m_arMinor[CoordinateSystem.utmk.rawValue] = 6356752.3141403558
        
        datum_params[0] = -146.43
        datum_params[1] = 507.89
        datum_params[2] = 681.46
        
        var tmp = m_arMinor[CoordinateSystem.wgs84.rawValue] / m_arMajor[CoordinateSystem.wgs84.rawValue]
        m_Es[CoordinateSystem.wgs84.rawValue] = 1.0 - tmp * tmp
        m_Esp[CoordinateSystem.wgs84.rawValue] = m_Es[CoordinateSystem.wgs84.rawValue] / (1.0 - m_Es[CoordinateSystem.wgs84.rawValue])
        
        if (m_Es[CoordinateSystem.wgs84.rawValue] < 0.00001) {
            m_Ind[CoordinateSystem.wgs84.rawValue] = 1.0
        } else {
            m_Ind[CoordinateSystem.wgs84.rawValue] = 0.0
        }
        
        tmp = m_arMinor[CoordinateSystem.katec.rawValue] / m_arMajor[CoordinateSystem.katec.rawValue]
        m_Es[CoordinateSystem.katec.rawValue] = 1.0 - tmp * tmp
        m_Esp[CoordinateSystem.katec.rawValue] = m_Es[CoordinateSystem.katec.rawValue] / (1.0 - m_Es[CoordinateSystem.katec.rawValue])
        
        if (m_Es[CoordinateSystem.katec.rawValue] < 0.00001) {
            m_Ind[CoordinateSystem.katec.rawValue] = 1.0
        } else {
            m_Ind[CoordinateSystem.katec.rawValue] = 0.0
        }
        
        tmp = m_arMinor[CoordinateSystem.tm.rawValue] / m_arMajor[CoordinateSystem.tm.rawValue]
        m_Es[CoordinateSystem.tm.rawValue] = 1.0 - tmp * tmp
        m_Esp[CoordinateSystem.tm.rawValue] = m_Es[CoordinateSystem.tm.rawValue] / (1.0 - m_Es[CoordinateSystem.tm.rawValue])
        
        if (m_Es[CoordinateSystem.tm.rawValue] < 0.00001) {
            m_Ind[CoordinateSystem.tm.rawValue] = 1.0
        } else {
            m_Ind[CoordinateSystem.tm.rawValue] = 0.0
        }
        
        tmp = m_arMinor[CoordinateSystem.utmk.rawValue] / m_arMajor[CoordinateSystem.utmk.rawValue]
        m_Es[CoordinateSystem.utmk.rawValue] = 1.0 - tmp * tmp
        m_Esp[CoordinateSystem.utmk.rawValue] = m_Es[CoordinateSystem.utmk.rawValue] / (1.0 - m_Es[CoordinateSystem.utmk.rawValue])
        
        if (m_Es[CoordinateSystem.utmk.rawValue] < 0.00001) {
            m_Ind[CoordinateSystem.utmk.rawValue] = 1.0
        } else {
            m_Ind[CoordinateSystem.utmk.rawValue] = 0.0
        }
        
        tmp = m_arMinor[CoordinateSystem.grs80.rawValue] / m_arMajor[CoordinateSystem.grs80.rawValue]
        m_Es[CoordinateSystem.grs80.rawValue] = 1.0 - tmp * tmp
        m_Esp[CoordinateSystem.grs80.rawValue] = m_Es[CoordinateSystem.grs80.rawValue] / (1.0 - m_Es[CoordinateSystem.grs80.rawValue])
        
        if (m_Es[CoordinateSystem.grs80.rawValue] < 0.00001) {
            m_Ind[CoordinateSystem.grs80.rawValue] = 1.0
        } else {
            m_Ind[CoordinateSystem.grs80.rawValue] = 0.0
        }
        
        src_m[CoordinateSystem.wgs84.rawValue] = m_arMajor[CoordinateSystem.wgs84.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.wgs84.rawValue]), e1fn(m_Es[CoordinateSystem.wgs84.rawValue]), e2fn(m_Es[CoordinateSystem.wgs84.rawValue]), e3fn(m_Es[CoordinateSystem.wgs84.rawValue]), m_arLatCenter[CoordinateSystem.wgs84.rawValue])
        dst_m[CoordinateSystem.wgs84.rawValue] = m_arMajor[CoordinateSystem.wgs84.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.wgs84.rawValue]), e1fn(m_Es[CoordinateSystem.wgs84.rawValue]), e2fn(m_Es[CoordinateSystem.wgs84.rawValue]), e3fn(m_Es[CoordinateSystem.wgs84.rawValue]), m_arLatCenter[CoordinateSystem.wgs84.rawValue])
        src_m[CoordinateSystem.katec.rawValue] = m_arMajor[CoordinateSystem.katec.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.katec.rawValue]), e1fn(m_Es[CoordinateSystem.katec.rawValue]), e2fn(m_Es[CoordinateSystem.katec.rawValue]), e3fn(m_Es[CoordinateSystem.katec.rawValue]), m_arLatCenter[CoordinateSystem.katec.rawValue])
        dst_m[CoordinateSystem.katec.rawValue] = m_arMajor[CoordinateSystem.katec.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.katec.rawValue]), e1fn(m_Es[CoordinateSystem.katec.rawValue]), e2fn(m_Es[CoordinateSystem.katec.rawValue]), e3fn(m_Es[CoordinateSystem.katec.rawValue]), m_arLatCenter[CoordinateSystem.katec.rawValue])
        src_m[CoordinateSystem.tm.rawValue] = m_arMajor[CoordinateSystem.tm.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.tm.rawValue]), e1fn(m_Es[CoordinateSystem.tm.rawValue]), e2fn(m_Es[CoordinateSystem.tm.rawValue]), e3fn(m_Es[CoordinateSystem.tm.rawValue]), m_arLatCenter[CoordinateSystem.tm.rawValue])
        dst_m[CoordinateSystem.tm.rawValue] = m_arMajor[CoordinateSystem.tm.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.tm.rawValue]), e1fn(m_Es[CoordinateSystem.tm.rawValue]), e2fn(m_Es[CoordinateSystem.tm.rawValue]), e3fn(m_Es[CoordinateSystem.tm.rawValue]), m_arLatCenter[CoordinateSystem.tm.rawValue])
        src_m[CoordinateSystem.grs80.rawValue] = m_arMajor[CoordinateSystem.grs80.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.grs80.rawValue]), e1fn(m_Es[CoordinateSystem.grs80.rawValue]), e2fn(m_Es[CoordinateSystem.grs80.rawValue]), e3fn(m_Es[CoordinateSystem.grs80.rawValue]), m_arLatCenter[CoordinateSystem.grs80.rawValue])
        dst_m[CoordinateSystem.grs80.rawValue] = m_arMajor[CoordinateSystem.grs80.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.grs80.rawValue]), e1fn(m_Es[CoordinateSystem.grs80.rawValue]), e2fn(m_Es[CoordinateSystem.grs80.rawValue]), e3fn(m_Es[CoordinateSystem.grs80.rawValue]), m_arLatCenter[CoordinateSystem.grs80.rawValue])
        src_m[CoordinateSystem.utmk.rawValue] = m_arMajor[CoordinateSystem.utmk.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.utmk.rawValue]), e1fn(m_Es[CoordinateSystem.utmk.rawValue]), e2fn(m_Es[CoordinateSystem.utmk.rawValue]), e3fn(m_Es[CoordinateSystem.utmk.rawValue]), m_arLatCenter[CoordinateSystem.utmk.rawValue])
        dst_m[CoordinateSystem.utmk.rawValue] = m_arMajor[CoordinateSystem.utmk.rawValue] * mlfn(e0fn(m_Es[CoordinateSystem.utmk.rawValue]), e1fn(m_Es[CoordinateSystem.utmk.rawValue]), e2fn(m_Es[CoordinateSystem.utmk.rawValue]), e3fn(m_Es[CoordinateSystem.utmk.rawValue]), m_arLatCenter[CoordinateSystem.utmk.rawValue])
    }
    
    private func D2R(_ degree: Double) -> Double {
        degree * .pi / 180.0
    }
    
    private func R2D(_ radian: Double) -> Double {
        radian * 180.0 / .pi
    }
    
    private func e0fn(_ x: Double) -> Double {
        1.0 - 0.25 * x * (1.0 + x / 16.0 * (3.0 + 1.25 * x))
    }
    
    private func e1fn(_ x: Double) -> Double {
        0.375 * x * (1.0 + 0.25 * x * (1.0 + 0.46875 * x))
    }
    
    private func e2fn(_ x: Double) -> Double {
        0.05859375 * x * x * (1.0 + 0.75 * x)
    }
    
    private func e3fn(_ x: Double) -> Double {
        return x * x * x * (35.0 / 3072.0)
    }
    
    private func mlfn(_ e0: Double, _ e1: Double, _ e2: Double, _ e3: Double, _ phi: Double) -> Double {
        e0 * phi - e1 * sin(2.0 * phi) + e2 * sin(4.0 * phi) - e3 * sin(6.0 * phi)
    }
    
    private func asinz(_ value: Double) -> Double {
        var value = value
        if (abs(value) > 1.0) { value = (value > 0 ? 1 : -1) }
        return asin(value)
    }
    
    func convert(from: CoordinateSystem, to: CoordinateSystem, geoPoint: GeoPoint) -> GeoPoint {
        let tmpPt = GeoPoint()
        let out_pt = GeoPoint()
        
        if (from == .wgs84) {
            tmpPt.longitude = D2R(geoPoint.longitude)
            tmpPt.latitude = D2R(geoPoint.latitude)
        } else {
            tm2geo(from, geoPoint, tmpPt)
        }
        
        if (to == .wgs84) {
            out_pt.longitude = R2D(tmpPt.longitude)
            out_pt.latitude = R2D(tmpPt.latitude)
        } else {
            geo2tm(to, tmpPt, out_pt)
            //out_pt.longitude = round(out_pt.longitude)
            //out_pt.latitude = round(out_pt.latitude)
        }
        
        return out_pt
    }
    
    private func geo2tm(_ dsttype: CoordinateSystem, _ in_pt: GeoPoint, _ out_pt: GeoPoint) {
        transform(.wgs84, dsttype, in_pt)
        let delta_lon = in_pt.longitude - m_arLonCenter[dsttype.rawValue]
        let sin_phi = sin(in_pt.latitude)
        let cos_phi = cos(in_pt.latitude)
        
        if (m_Ind[dsttype.rawValue] != 0) {
            let b = cos_phi * sin(delta_lon)
            
            if ((abs(abs(b) - 1.0)) < EPSLN) {
                //Log.d("무한대 에러")
                //System.out.println("무한대 에러")
            }
        } else {
            let b = 0.0
            var con = acos(cos_phi * cos(delta_lon) / sqrt(1.0 - b * b))
            
            if (in_pt.latitude < 0) {
                con = con * -1
            }
        }
        
        let al = cos_phi * delta_lon
        let als = al * al
        let c = m_Esp[dsttype.rawValue] * cos_phi * cos_phi
        let tq = tan(in_pt.latitude)
        let t = tq * tq
        let con = 1.0 - m_Es[dsttype.rawValue] * sin_phi * sin_phi
        let n = m_arMajor[dsttype.rawValue] / sqrt(con)
        let ml = m_arMajor[dsttype.rawValue] * mlfn(e0fn(m_Es[dsttype.rawValue]), e1fn(m_Es[dsttype.rawValue]), e2fn(m_Es[dsttype.rawValue]), e3fn(m_Es[dsttype.rawValue]), in_pt.latitude)
        
        out_pt.longitude = m_arScaleFactor[dsttype.rawValue] * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 * (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * m_Esp[dsttype.rawValue]))) + m_arFalseEasting[dsttype.rawValue]
        out_pt.latitude = m_arScaleFactor[dsttype.rawValue] * (ml - dst_m[dsttype.rawValue] + n * tq * (als * (0.5 + als / 24.0 * (5.0 - t + 9.0 * c + 4.0 * c * c + als / 30.0 * (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * m_Esp[dsttype.rawValue]))))) + m_arFalseNorthing[dsttype.rawValue]
    }
    
    private func tm2geo(_ srctype: CoordinateSystem, _ in_pt: GeoPoint, _ out_pt: GeoPoint) {
        let tmpPt = GeoPoint(longitude: in_pt.longitude, latitude: in_pt.latitude)
        let max_iter = 6
        
        if (m_Ind[srctype.rawValue] != 0) {
            let f = exp(in_pt.longitude / (m_arMajor[srctype.rawValue] * m_arScaleFactor[srctype.rawValue]))
            let g = 0.5 * (f - 1.0 / f)
            let temp = m_arLatCenter[srctype.rawValue] + tmpPt.latitude / (m_arMajor[srctype.rawValue] * m_arScaleFactor[srctype.rawValue])
            let h = cos(temp)
            let con = sqrt((1.0 - h * h) / (1.0 + g * g))
            out_pt.latitude = asinz(con)
            
            if (temp < 0) { out_pt.latitude *= -1 }
            
            if ((g == 0) && (h == 0)) {
                out_pt.longitude = m_arLonCenter[srctype.rawValue]
            } else {
                out_pt.longitude = atan(g / h) + m_arLonCenter[srctype.rawValue]
            }
        }
        
        tmpPt.longitude -= m_arFalseEasting[srctype.rawValue]
        tmpPt.latitude -= m_arFalseNorthing[srctype.rawValue]
        
        let con = (src_m[srctype.rawValue] + tmpPt.latitude / m_arScaleFactor[srctype.rawValue]) / m_arMajor[srctype.rawValue]
        var phi = con
        
        var i = 0
        
        while (true) {
            let delta_Phi = ((con + e1fn(m_Es[srctype.rawValue]) * sin(2.0 * phi) - e2fn(m_Es[srctype.rawValue]) * sin(4.0 * phi) + e3fn(m_Es[srctype.rawValue]) * sin(6.0 * phi)) / e0fn(m_Es[srctype.rawValue])) - phi
            phi = phi + delta_Phi
            
            if (abs(delta_Phi) <= EPSLN) { break }
            
            if (i >= max_iter) {
                //Log.d("무한대 에러")
                //System.out.println("무한대 에러")
                break
            }
            
            i += 1
        }
        
        if (abs(phi) < (.pi / 2)) {
            let sin_phi = sin(phi)
            let cos_phi = cos(phi)
            let tan_phi = tan(phi)
            let c = m_Esp[srctype.rawValue] * cos_phi * cos_phi
            let cs = c * c
            let t = tan_phi * tan_phi
            let ts = t * t
            let cont = 1.0 - m_Es[srctype.rawValue] * sin_phi * sin_phi
            let n = m_arMajor[srctype.rawValue] / sqrt(cont)
            let r = n * (1.0 - m_Es[srctype.rawValue]) / cont
            let d = tmpPt.longitude / (n * m_arScaleFactor[srctype.rawValue])
            let ds = d * d
            out_pt.latitude = phi - (n * tan_phi * ds / r) * (0.5 - ds / 24.0 * (5.0 + 3.0 * t + 10.0 * c - 4.0 * cs - 9.0 * m_Esp[srctype.rawValue] - ds / 30.0 * (61.0 + 90.0 * t + 298.0 * c + 45.0 * ts - 252.0 * m_Esp[srctype.rawValue] - 3.0 * cs)))
            out_pt.longitude = m_arLonCenter[srctype.rawValue] + (d * (1.0 - ds / 6.0 * (1.0 + 2.0 * t + c - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * m_Esp[srctype.rawValue] + 24.0 * ts))) / cos_phi)
        } else {
            out_pt.latitude = .pi * 0.5 * sin(tmpPt.latitude)
            out_pt.longitude = m_arLonCenter[srctype.rawValue]
        }
        transform(srctype, .wgs84, out_pt)
    }
    
    func getDistancebyGeo(_ pt1: GeoPoint , _ pt2: GeoPoint) -> Double {
        let lat1 = D2R(pt1.latitude)
        let lon1 = D2R(pt1.longitude)
        let lat2 = D2R(pt2.latitude)
        let lon2 = D2R(pt2.longitude)
        
        let longitude = lon2 - lon1
        let latitude = lat2 - lat1
        
        let a = pow(sin(latitude / 2.0), 2) + cos(lat1) * cos(lat2) * pow(sin(longitude / 2.0), 2)
        return 6376.5 * 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    }
    
    /**
     * convert between geodetic coordinates (longitude, latitude, height)
     * and gecentric coordinates (X, Y, Z)
     * ported from Proj 4.9.9 geocent.c
     */
    
    private let HALF_PI = 0.5 * .pi
    private let COS_67P5 = 0.38268343236508977  /* cosine of 67.5 degrees */
    private let AD_C = 1.0026000
    
    private func transform(_ srctype: CoordinateSystem, _ dsttype: CoordinateSystem, _ point: GeoPoint) {
        if (srctype == dsttype) {
            return
        }
        
        if ((srctype != .wgs84 && srctype != .grs80 && srctype != .utmk) || (dsttype != .wgs84 && dsttype != .grs80 && dsttype != .utmk)) {
            // Convert to geocentric coordinates.
            geodetic_to_geocentric(srctype, point)
            
            // Convert between datums
            if (srctype != .wgs84 && srctype != .grs80 && srctype != .utmk) {
                geocentric_to_wgs84(point)
            }
            
            if (dsttype != .wgs84 && dsttype != .grs80 && dsttype != .utmk) {
                geocentric_from_wgs84(point)
            }
            
            // Convert back to geodetic coordinates
            geocentric_to_geodetic(dsttype, point)
        }
    }
    
    private func geodetic_to_geocentric(_ type: CoordinateSystem, _ p: GeoPoint) -> Bool {
        /*
         * The function Convert_Geodetic_To_Geocentric converts geodetic coordinates
         * (latitude, longitude, and height) to geocentric coordinates (X, Y, Z),
         * according to the current ellipsoid parameters.
         *
         *    Latitude  : Geodetic latitude in radians                     (input)
         *    Longitude : Geodetic longitude in radians                    (input)
         *    Height    : Geodetic height, in meters                       (input)
         *    X         : Calculated Geocentric X coordinate, in meters    (output)
         *    Y         : Calculated Geocentric Y coordinate, in meters    (output)
         *    Z         : Calculated Geocentric Z coordinate, in meters    (output)
         *
         */
        
        var Longitude = p.longitude
        var Latitude = p.latitude
        let Height = p.altitude
        var X = 0.0  // output
        var Y = 0.0
        var Z = 0.0
        
        var Rn = 0.0            /*  Earth radius at location  */
        var Sin_Lat = 0.0       /*  sin(Latitude)  */
        var Sin2_Lat = 0.0      /*  Square of sin(Latitude)  */
        var Cos_Lat = 0.0       /*  cos(Latitude)  */
        
        /*
         ** Don't blow up if Latitude is just a little out of the value
         ** range as it may just be a rounding issue.  Also removed longitude
         ** test, it should be wrapped by cos() and sin().  NFW for PROJ.4, Sep/2001.
         */
        if (Latitude < -HALF_PI && Latitude > -1.001 * HALF_PI ) {
            Latitude = -HALF_PI
        } else if (Latitude > HALF_PI && Latitude < 1.001 * HALF_PI ) {
            Latitude = HALF_PI
        } else if ((Latitude < -HALF_PI) || (Latitude > HALF_PI)) { /* Latitude out of range */
            return true
        }
        
        /* no errors */
        if (Longitude > .pi) {
            Longitude -= (2 * .pi)
        }
        Sin_Lat = sin(Latitude)
        Cos_Lat = cos(Latitude)
        Sin2_Lat = Sin_Lat * Sin_Lat
        Rn = m_arMajor[type.rawValue] / (sqrt(1.0e0 - m_Es[type.rawValue] * Sin2_Lat))
        X = (Rn + Height) * Cos_Lat * cos(Longitude)
        Y = (Rn + Height) * Cos_Lat * sin(Longitude)
        Z = ((Rn * (1 - m_Es[type.rawValue])) + Height) * Sin_Lat
        
        p.longitude = X
        p.latitude = Y
        p.altitude = Z
        return false
    } // cs_geodetic_to_geocentric()
    
    /** Convert_Geocentric_To_Geodetic
     * The method used here is derived from 'An Improved Algorithm for
     * Geocentric to Geodetic Coordinate Conversion', by Ralph Toms, Feb 1996
     */
    private func geocentric_to_geodetic(_ type: CoordinateSystem, _ p: GeoPoint) {
        let X = p.longitude
        let Y = p.latitude
        let Z = p.altitude
        var Longitude = 0.0
        var Latitude = 0.0
        var Height = 0.0
        
        var W = 0.0        /* distance from Z axis */
        var W2 = 0.0       /* square of distance from Z axis */
        var T0 = 0.0       /* initial estimate of vertical component */
        var T1 = 0.0       /* corrected estimate of vertical component */
        var S0 = 0.0       /* initial estimate of horizontal component */
        var S1 = 0.0       /* corrected estimate of horizontal component */
        var Sin_B0 = 0.0   /* sin(B0), B0 is estimate of Bowring aux doubleiable */
        var Sin3_B0 = 0.0  /* cube of sin(B0) */
        var Cos_B0 = 0.0   /* cos(B0) */
        var Sin_p1 = 0.0   /* sin(phi1), phi1 is estimated latitude */
        var Cos_p1 = 0.0   /* cos(phi1) */
        var Rn = 0.0       /* Earth radius at location */
        var Sum = 0.0      /* numerator of cos(phi1) */
        var At_Pole = false  /* indicates location is in polar region */
        
        if (X != 0.0) {
            Longitude = atan2(Y,X)
        }
        else {
            if (Y > 0) {
                Longitude = HALF_PI
            }
            else if (Y < 0) {
                Longitude = -HALF_PI
            }
            else {
                At_Pole = true
                Longitude = 0.0
                if (Z > 0.0) {  /* north pole */
                    Latitude = HALF_PI
                }
                else if (Z < 0.0) {  /* south pole */
                    Latitude = -HALF_PI
                }
                else {  /* center of earth */
                    Latitude = HALF_PI
                    Height = -m_arMinor[type.rawValue]
                    return
                }
            }
        }
        W2 = X*X + Y*Y
        W = sqrt(W2)
        T0 = Z * AD_C
        S0 = sqrt(T0 * T0 + W2)
        Sin_B0 = T0 / S0
        Cos_B0 = W / S0
        Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0
        T1 = Z + m_arMinor[type.rawValue] * m_Esp[type.rawValue] * Sin3_B0
        Sum = W - m_arMajor[type.rawValue] * m_Es[type.rawValue] * Cos_B0 * Cos_B0 * Cos_B0
        S1 = sqrt(T1*T1 + Sum * Sum)
        Sin_p1 = T1 / S1
        Cos_p1 = Sum / S1
        Rn = m_arMajor[type.rawValue] / sqrt(1.0 - m_Es[type.rawValue] * Sin_p1 * Sin_p1)
        if (Cos_p1 >= COS_67P5) {
            Height = W / Cos_p1 - Rn
        }
        else if (Cos_p1 <= -COS_67P5) {
            Height = W / -Cos_p1 - Rn
        }
        else {
            Height = Z / Sin_p1 + Rn * (m_Es[type.rawValue] - 1.0)
        }
        if (At_Pole == false) {
            Latitude = atan(Sin_p1 / Cos_p1)
        }
        
        p.longitude = Longitude
        p.latitude = Latitude
        p.altitude = Height
    } // geocentric_to_geodetic()
    
    /****************************************************************/
    // geocentic_to_wgs84(defn, p )
    //  defn = coordinate system definition,
    //  p = point to transform in geocentric coordinates (x,y,z)
    private func geocentric_to_wgs84(_ p: GeoPoint) {
        p.longitude += datum_params[0]
        p.latitude += datum_params[1]
        p.altitude += datum_params[2]
    } // geocentric_to_wgs84
    
    /****************************************************************/
    // geocentic_from_wgs84()
    //  coordinate system definition,
    //  point to transform in geocentric coordinates (x,y,z)
    private func geocentric_from_wgs84(_ p: GeoPoint) {
        p.longitude -= datum_params[0]
        p.latitude -= datum_params[1]
        p.altitude -= datum_params[2]
    } //geocentric_from_wgs84()
}
